{-=FastaRegionRandomizer (FRR): A Haskell-based solution to=-}
{-=generate X random snvs given chromosome, sequence window string,=-}
{-=strand, number of runs per batch and a fasta file.=-}
{-=Author: Matthew Mosior=-}
{-=Version: 2.0=-}
{-=Synopsis:  This Haskell Script will take in=-}
{-=user specified chromsome, start, stop, strand,=-} 
{-=number of runs per batch, and fasta file=-}
{-=and will generate X random SNVs.=-}


{-Lanuguage Extension.-}

{-# LANGUAGE MultiWayIf #-}

{----------------------}


{-Imports-}

import Bio.Core.Sequence as BCS
import Bio.Sequence.Fasta as BSF
import Control.Arrow as CA
import Control.Monad as CM
import Control.Monad.Primitive as CMP
import Data.ByteString as DB
import Data.ByteString.Char8 as DBC
import Data.ByteString.Lazy as DBL
import Data.ByteString.Lazy.Char8 as DBLC8
import Data.Char as DC
import Data.Foldable as DF
import Data.Int as DI
import Data.List as DL
import Data.List.Split as DLS
import Data.Ord as DO
import Data.Traversable as DT
import Data.Vector.Fusion.Stream.Monadic as DVFSM
import Data.Word8 as DW8
import System.Console.GetOpt as SCG
import System.Process as SP
import System.Random as SR
import System.Random.MWC as SRMWC
import System.Environment as SE
import System.Exit as SX
import System.IO as SIO
import System.IO.Temp as SIOT

{---------}


{-Custom CML Option Datatype.-}

data Flag
    = Verbose           -- -v
    | Version           -- -V -?
    | BatchSize  String -- -b
    | OutputFile String -- -o
    | Help              -- --help
    deriving (Eq,Ord,Show)

{-----------------------------}


{-Custom bool functions for Flag Datatype.-}

--isOutputFile -> This function will
--test for OutputFile flag.
isOutputFile :: Flag -> Bool
isOutputFile (OutputFile _) = True
isOutputFile _              = False

--isBatchSize -> This function will
--test for BatchSize flag.
isBatchSize :: Flag -> Bool
isBatchSize (BatchSize _) = True
isBatchSize _             = False

{------------------------------------------}


{-Custom extraction functions for Flag Datatype.-}

--extractOutputFile -> This function will
--extract the string associated with 
--OutputFile.
extractOutputFile :: Flag -> String
extractOutputFile (OutputFile x) = x

--extractBatchSize -> This function will
--extract the string associated with
--BatchSize.
extractBatchSize :: Flag -> String
extractBatchSize (BatchSize x) = x

{------------------------------------------------}


{-Custom extraction functions for Sequence Datatype.-}

--extractSeqData -> This function will
--extract the SeqData field of Seq constructor 
--of the Sequence Datatype.
extractSeqData :: Sequence -> SeqData
extractSeqData (Seq _ x _) = x

--extractSeqLabel -> This function will
--extract the SeqLabel field of Seq constructor
--of the Sequence Datatype.
extractSeqLabel :: Sequence -> SeqLabel
extractSeqLabel (Seq x _ _) = x

--extractQualData -> This function will
--extract the QualData field of Seq constructor
--of the Sequence Datatype.
extractQualData :: Sequence -> Maybe QualData
extractQualData (Seq _ _ x) = x

--extractunSD -> This function will
--extract the unSD field of SeqData constructor.
extractunSD :: SeqData -> DBL.ByteString
extractunSD (SeqData unSD) = unSD

--extractunSL -> This function will
--extract the unSL field of SeqLabel constructor.
extractunSL :: SeqLabel -> DBL.ByteString
extractunSL (SeqLabel unSL) = unSL

--extractunQD -> This function will
--extract the unQD field of QualData constructor.
extractunQD :: QualData -> DBL.ByteString
extractunQD (QualData unQD) = unQD

{----------------------------------------------------}


{-Option Description function relating to datatype above.-}

--options -> This function will
--describe flags.
options :: [OptDescr Flag]
options =
    [ Option ['v']     ["verbose"]    (NoArg Verbose)                "Output on stderr.",
      Option ['V','?'] ["version"]    (NoArg Version)                "Show version number.",
      Option ['b']     ["batchsize"]  (ReqArg BatchSize "BATCHSIZE") "The number of batches (DEFAULT=1).",
      Option ['o']     ["outputfile"] (ReqArg OutputFile "OUTFILE")  "The path to the output file.", 
      Option []        ["help"]       (NoArg Help)                   "Print this help message."
    ]

--compilerOpts -> This function will
--parse incoming command line arguments.
compilerOpts :: [String] -> IO ([Flag],[String])
compilerOpts argv =
    case getOpt Permute Main.options argv of
        (args,files,[]) ->
            if | DL.elem Help args ->
               do SIO.hPutStrLn stderr (greeting DL.++ SCG.usageInfo header Main.options)
                  SX.exitWith SX.ExitSuccess
               | DL.elem Version args ->
               do SIO.hPutStrLn stderr (greeting DL.++ version DL.++ SCG.usageInfo header Main.options)
                  SX.exitWith SX.ExitSuccess
               | (DL.length (DL.filter (isBatchSize) args) > 0) &&
                 (not (intCheck (extractBatchSize (DL.head (DL.filter (isBatchSize) args))))) ->
               do SIO.hPutStrLn stderr (bserror DL.++ github DL.++ SCG.usageInfo header Main.options)
                  SX.exitWith (SX.ExitFailure 1)
               | (DL.length files > 3 || DL.length files < 3) ->
               do SIO.hPutStrLn stderr (flerror DL.++ github DL.++ SCG.usageInfo header Main.options)
                  SX.exitWith (SX.ExitFailure 1)
               | otherwise ->
               return (DL.nub args, files)
        (_,_,errors) -> do
            SIO.hPutStrLn stderr (DL.concat errors DL.++ SCG.usageInfo header Main.options)
            SX.exitWith (SX.ExitFailure 1)
        where
            greeting        = "\nFasta Region Randomizer, Copyright (c) 2020 Matthew Mosior.\n"
            header          = "\nUsage: frr [-vV?bo] [Sequence Window String] [Total Number of Randomized Variants] [Fasta file]"
            version         = "\nFasta Region Randomizer (FRR), Version 1.0.\n"
            github          = "\nPlease see https://github.com/Matthew-Mosior/Fasta-Region-Randomizer/wiki for more information.\n"
            flerror         = "\nIncorrect number of input arguments:\n\
                               \Please provide exactly three inputs:\n\
                               \First argument  -> Sequence Window String\n\ 
                               \Second argument -> Total number of randomized variants to be created\n\
                               \Third argument  -> Fasta file\n"
            bserror         = "\nIncorrect batch size argument.\n\
                               \Please provide a number (integer).\n"

{---------------------------------------------------------}

{-General Utility Functions.-}

--tuplifyTwo -> The function will
--turn a list of length two to a tuple.
tuplifyTwo :: [a] -> (a,a)
tuplifyTwo [x,y] = (x,y)

--tripletThrd -> This function will
--grab the third element of a triplet.
tripletThrd :: (a,a,a) -> a
tripletThrd (_,_,a) = a

--mapNotLast -> This function will
--work like the traditional map 
--function in Data.List, but not
--map to the last element of a list.
mapNotLast :: (a -> a) -> [a] -> [a]
mapNotLast fn []     = []
mapNotLast fn [x]    = [x]
mapNotLast fn (x:xs) = fn x : mapNotLast fn xs

--bslToStr -> This function will
--Convert from Bytestring (Lazy) to String.
bslToStr :: DBL.ByteString -> String
bslToStr = DL.map (DC.chr . fromEnum) . DBL.unpack

--bslToStrSingleton -> This function will
--convert from a singleton Bytestring (Lazy) to String.
bslToStrS :: Word8 -> String
bslToStrS = DBLC8.unpack . DBL.singleton

--strToBSC8 -> This function will
--convert Strings to Bytestring (Char8).
strToBSC8 :: String -> DBC.ByteString
strToBSC8 xs = DBC.pack xs

--atRandomIndex -> This function will
--return random element of a list.
atRandomIndex :: [a] -> IO a
atRandomIndex xs = do
    index <- SR.randomRIO (0,DL.length xs-1)
    return $ xs DL.!! index

{----------------------------}


{-Batch Size functions.-}

--intCheck -> This function will
--check that the batch size string
--provided by the user is made up
--a number (integer).
intCheck :: String -> Bool
intCheck = DL.all DC.isDigit

{-----------------------}


{-Random positions function.-}

--allRandomIndices -> This function will
--generate random indices from 1 to the
--length of the sequence window string
--list X times, where X is the total number
--of randomized variants.
allRandomIndices :: [(String,Int,Char)] -> String -> IO [(String,Int,Char)]
allRandomIndices _  [] = return []
allRandomIndices [] _  = return []
allRandomIndices xs ys = do
    --Create the random indices.
    let randomindicesmonadicstream = DVFSM.replicateM (Prelude.read ys) (atRandomIndex xs)
    --Convert the monadic stream to a list.
    DVFSM.toList randomindicesmonadicstream

{----------------------------}


{-Custom Sequence Data Structure Functions.-}

--This function will
--create a custom sequence based 
--on the sequence window string.
--[(String,Int,Char)] <-> [(DescriptionLine,Position,nucleotide/aminoacid sequence)]
createCustomSequence :: [String] -> [Sequence] -> [(String,Int,Char)]
createCustomSequence _      [] = []
createCustomSequence []     _  = []
createCustomSequence (x:xs) ys = zipped DL.++ (createCustomSequence xs ys)  
    where
        --Local definitions.--
        zipped = DL.zip3 alldlines allpositions sequencechunk
        allpositions = [start..stop]
        start = read (fst (tuplifyTwo (DLS.splitOn "-" (DL.concat (DL.tail (DLS.splitOn ":" x))))))
        stop  = read (snd (tuplifyTwo (DLS.splitOn "-" (DL.concat (DL.tail (DLS.splitOn ":" x))))))
        alldlines = DL.replicate (DL.length (DBLC8.unpack
                                            (grabChunkByIndices (grabFastaSequence
                                                                (DL.head (DLS.splitOn ":" x)) ys)
                                                                (fromIntegral start :: Int64
                                                                ,fromIntegral stop :: Int64))))
                                 (DL.head (DLS.splitOn ":" x))
        sequencechunk = DBLC8.unpack (grabChunkByIndices (grabFastaSequence
                                                         (DL.head (DLS.splitOn ":" x)) ys)
                                                         (fromIntegral start :: Int64
                                                         ,fromIntegral stop :: Int64))
        ----------------------

--grabFastaSequence -> This function will
--grab the correct fasta sequence
--using chromosome information
--in the region file.
grabFastaSequence :: String -> [Sequence] -> DBL.ByteString
grabFastaSequence x ys = smallGrabFastaSequence x ys [0..(DL.length ys) - 1]

--smallGrabFastaSequence -> This function will
--grab the correct fasta sequence
--using chromosome information
--in the region file.
smallGrabFastaSequence :: String -> [Sequence] -> [Int] -> DBL.ByteString
smallGrabFastaSequence _ _ [] = DBL.empty
smallGrabFastaSequence x ys (z:zs) = if | ((bslToStr (extractunSL (extractSeqLabel (ys DL.!! z)))) == x) ->
                                        extractunSD (extractSeqData (ys DL.!! z))
                                        | otherwise ->
                                        smallGrabFastaSequence x ys zs

--grabChunkByIndices -> This function will
--grab a chunk of a ByteString Based on
--start and end indices.
grabChunkByIndices :: DBLC8.ByteString -> (Int64,Int64) -> DBLC8.ByteString
grabChunkByIndices xs (a,b) = DBLC8.take ((b-1) - (a-1)) (DBLC8.drop (a-1) xs)

{--------------------------------------------------}


{-Random nucleotides function.-}

--randomNucleotides -> This function will 
--create a list of random integers from 0 
--to 2.
randomNucleotides :: GenIO -> String -> IO [Int]
randomNucleotides _   [] = return []
randomNucleotides gen xs = do
    --Create the random nucleotides using a position mapping.
    let randomnucleotidesmonadicstream = DVFSM.replicateM (Prelude.read xs) ((SRMWC.uniformR (0,2) gen) :: IO Int)
    --Convert randomnucleotidesmonadicstream to a list.
    DVFSM.toList randomnucleotidesmonadicstream

{------------------------------}


{-Random snv generator functions.-}

--randomMapping -> This function will
--map a nucleotide to a nucleotide
--transistion/transversion.
randomMapping :: [(String,Int,Char)] -> [Int] -> [[String]]
randomMapping [] _  = []
randomMapping _  [] = []
randomMapping xs ys = DL.map (\(a,b,c,(d,e)) -> [a,b,c,d,e]) (randomMappingSmall xs ys)

--randomMappingSmall -> This function will
--map a nucleotide to a nucleotide
--transistion/transversion.
randomMappingSmall :: [(String,Int,Char)] -> [Int] -> [(String,String,String,(String,String))]
randomMappingSmall []     _  = []
randomMappingSmall _      [] = []
randomMappingSmall (x:xs) (y:ys) = [((\(a,_,_) -> a) x
                                    ,Prelude.show ((\(_,b,_) -> b) x)
                                    ,Prelude.show ((\(_,b,_) -> b) x)
                                    ,([(\(_,_,c) -> c) x]
                                    ,[snd ((DL.concatMap (\s -> DL.filter 
                                                         (\(r,_) -> r == ((\(_,_,c) -> c) s)) 
                                                         nucleotidemapping) [x]) DL.!! y)]))] 
                                    DL.++ (randomMappingSmall xs ys) 

--nucleotidemapping -> List containing possible mappings for nucleotides mutations.
--(NUCLEOTIDE,NUCLEOTIDE_TRANSITION/NUCLEOTIDE_TRANSVERSION)
nucleotidemapping :: [(Char,Char)]
nucleotidemapping = [('A','T'),('A','G'),('A','C')
                    ,('T','A'),('T','G'),('T','C')
                    ,('G','C'),('G','A'),('G','T')
                    ,('C','G'),('C','A'),('C','T')]

{-----------------------------------}


{-Batch adder function.-}

--batchAdder -> This function will
--add column to each list
--annotating the batch.
batchAdder :: [[String]] -> String -> [Flag] -> [[String]]
batchAdder [] [] []   = []
batchAdder xs ys opts = --Check if user passed BatchSize flag.
                        if DL.length (DL.filter (isBatchSize) opts) > 0
                            then DL.map (DL.concat) (DL.map (\(a,b) -> [a,b]) (DL.zip batchlist xs)) 
                            else DL.map (\x -> "1" : x) xs
    where
        --Local definitions.--
        --batchlist -> This list is what will be added to xs (used in then statement).
        batchlist = DL.map (\x -> [x]) 
                           (DL.map (Prelude.show) 
                           (DL.concat 
                           (DL.map 
                           (\x -> DL.take ((Prelude.read ys) `div` 
                           ((Prelude.read (extractBatchSize (DL.head (DL.filter (isBatchSize) opts))) :: Int))) 
                           (DL.repeat x)) 
                           [1..(Prelude.read (extractBatchSize (DL.head (DL.filter (isBatchSize) opts))) :: Int)])))
        ----------------------

{-----------------------}


{-Printing functions.-}

--tempFileCreation -> This function will
--print the file to stdout using
--readProcess of the unix tool cat.
catFile :: [[String]] -> IO ()
catFile [] = return ()
catFile xs = do
    --Open a temporary file.
    (tempfile,temph) <- SIOT.openTempFile "." "temp.txt"
    --Intercalate a tab, and then a newline into xs.
    let intercalatedxs = DL.intercalate "\n" (DL.map (DL.intercalate "\t") xs)
    --Add intercalatedxs to temp.txt.
    SIO.hPutStrLn temph intercalatedxs
    --Close the temporary file's handle.
    hClose temph
    --Print out the contents of tempfile to the screen using cat unix tool.
    (_,_,_,ph) <- SP.createProcess (SP.proc "cat" [tempfile])
    ec <- SP.waitForProcess ph
    case ec of
        SX.ExitSuccess   -> do _ <- SP.readProcess "rm" [tempfile] []
                               return ()
        SX.ExitFailure _ -> do _ <- error "Could not cat file."
                               _ <- SP.readProcess "rm" [tempfile] []
                               return ()

--printFile -> This function will
--print the file to either stdout
--or to a output file based on
--command-lines options provided.
printFile :: [Flag] -> [[String]] -> IO ()
printFile [] [] = return ()
printFile [] _  = return ()
printFile _  [] = return ()
printFile opts xs = do
    --Grab just "OUTFILE".
    let outfile = DL.head (DL.filter (isOutputFile) opts)
    --Extract the string from FilterFields.
    let outfilestring = extractOutputFile outfile
    --mapNotLast tabs and newlines in xs.
    let tabsandnewlinesadded = DL.intercalate "\n" (DL.map (DL.intercalate "\t") xs)
    --Write the output to the user-specified filename.
    SIO.writeFile (outfilestring) $ (tabsandnewlinesadded)

{---------------------}


{-FRR Specific Function.-}

--processArgsAndFiles -> This function will
--walk through each of the command-line
--arguments and files provided by the user.
processArgsAndFiles :: ([Flag],[String]) -> IO ()
processArgsAndFiles ([],[]) = return () 
processArgsAndFiles (options,files) = do
    --Process the Sequence Window String argument.
    let readseqwindowstrarg = files DL.!! 0
    ----------------------------------------------
    --Process the total number of randomized variants.
    let readtotalnumberofrandomvararg = files DL.!! 1
    --------------------------------------------------
    --Process the fasta file.
    readfastafile <- BSF.readFasta (files DL.!! 2)   
    -------------------------
    --Prepare readseqwindowstrarg.
    let seqwindowstr = DLS.splitOn "#" (DL.init (DL.tail readseqwindowstrarg)) 
    ------------------------------
    --Create custom sequence using sequence window string
    --and fasta file.
    let customsequence = createCustomSequence seqwindowstr readfastafile 
    --Generate all list of randomly picked sequence window
    --strings [Total Amount of Randomized Variants] times.
    randomindices <- allRandomIndices customsequence readtotalnumberofrandomvararg
    ------------------------------
    --Initialize random number generator by seeding the MWC PRNG.
    randgen <- SRMWC.createSystemRandom
    -------------------------------------------------------------
    --Get random nucleotides.
    randomnucleotides <- randomNucleotides randgen readtotalnumberofrandomvararg
    -------------------------
    --Map the current sequence to the randomized sequences
    --using a nucleotidemapping.
    let randomsnvs = randomMapping randomindices randomnucleotides
    ----------------------------
    --Add batch column to randomsnvs.
    let batchadded = batchAdder randomsnvs readtotalnumberofrandomvararg options 
    ---------------------------------
    --Add header to randomsnvs.
    let finalpreprint = [["Batch","Chromosome","Start","Stop","Ref","Alt"]] DL.++ batchadded
    ---------------------------
    --Print the file to stdout (cat) or to a file.
    if DL.length (DL.filter (isOutputFile) options) > 0 
        then printFile options finalpreprint
        else catFile finalpreprint
    ----------------------------------------------

{-------------------------}


{-Main function.-}

main :: IO ()
main = do
    --Get command line arguments.
    (args,files) <- SE.getArgs >>= compilerOpts
    --Run args and files through processArgsandFiles.
    processArgsAndFiles (args,files)

{----------------}
