{-=FastaRegionRandomizer (FRR): A Haskell-based solution to=-}
{-=generate X random snvs given chromosome, start, stop,=-}
{-=strand, number of runs per batch and a fasta file.=-}
{-=Author: Matthew Mosior=-}
{-=Version: 1.0=-}
{-=Synopsis:  This Haskell Script will take in=-}
{-=user specified chromsome, start, stop, strand,=-} 
{-=number of runs per batch, and fasta file=-}
{-=and will generate X random SNVs.=-}


{-Imports-}

import Bio.Core.Sequence as BCS
import Bio.Sequence.Fasta as BSF
import Control.Arrow as CA
import Control.Monad as CM
import Data.ByteString as DB
import Data.ByteString.Char8 as DBC
import Data.ByteString.Lazy as DBL
import Data.ByteString.Search.DFA as DBSDFA
import Data.Char as DC
import Data.List as DL
import Data.List.Split as DLS
import Data.Ord as DO
import Data.Traversable as DT
import System.Console.GetOpt as SCG
import System.Process as SP
import System.Random as SR
import System.Environment as SE
import System.Exit as SX
import System.IO as SIO
import System.IO.Temp as SIOT
import Text.PrettyPrint.Boxes as TPB

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
            if DL.elem Help args
                then do SIO.hPutStrLn stderr (greeting ++ SCG.usageInfo header Main.options)
                        SX.exitWith SX.ExitSuccess
                else if DL.elem Version args
                    then do SIO.hPutStrLn stderr (greeting ++ version ++ SCG.usageInfo header Main.options)
                            SX.exitWith SX.ExitSuccess
                    else if (DL.length (DL.filter (isBatchSize) args) > 0) &&
                            (not (intCheck (extractBatchSize (DL.head (DL.filter (isBatchSize) args)))))
                        then do SIO.hPutStrLn stderr (bserror ++ github ++ SCG.usageInfo header Main.options)
                                SX.exitWith (SX.ExitFailure 1)
                        else if (DL.length files > 5 || DL.length files < 5)
                            then do SIO.hPutStrLn stderr (flerror ++ github ++ SCG.usageInfo header Main.options)
                                    SX.exitWith (SX.ExitFailure 1)
                            else return (DL.nub args, files)
        (_,_,errors) -> do
            SIO.hPutStrLn stderr (DL.concat errors ++ SCG.usageInfo header Main.options)
            SX.exitWith (SX.ExitFailure 1)
        where
            greeting        = "Fasta Region Randomizer, Copyright (c) 2020 Matthew Mosior.\n"
            header          = "Usage: frr [-vV?o] [Chromosome] [Start] [Stop] [Number of runs per batch] [Fasta file]"
            version         = "Fasta Region Randomizer (FRR), Version 1.0.\n"
            github          = "Please see https://github.com/Matthew-Mosior/Fasta-Region-Randomizer/wiki for more information.\n"
            flerror         = "Incorrect number of input arguments:\n\
                              \Please provide exactly five arguments files.\n\
                              \First argument   -> Chromosome\n\
                              \Second argument  -> Start\n\
                              \Third argument   -> Stop\n\ 
                              \Fourth argument  -> Number of runs per batch\n\
                              \Fifth argument   -> Fasta file\n"
            bserror         = "Incorrect batch size argument.\n\
                              \Please provide a number (integer).\n"

{---------------------------------------------------------}

{-General Utility Functions.-}

--mapTuple -> This function will
--map a function across all elements
--of a two-tuple.
mapTuple = CM.join (***)

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

--strToBSC8 -> This function will
--convert Strings to Bytestring (Char8).
strToBSC8 :: String -> DBC.ByteString
strToBSC8 xs = DBC.pack xs

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

--randomPositions -> This function will
--generate random positions within
--specified start and stop.
randomPositions :: String -> String -> String -> [Flag] -> IO [Int]
randomPositions [] [] [] []   = return []
randomPositions _  [] [] []   = return []
randomPositions [] _  [] []   = return []
randomPositions [] [] _  []   = return []
randomPositions [] [] [] _    = return []
randomPositions xs ys zs opts = --Check if user passed BatchSize flag.
                                if DL.length (DL.filter (isBatchSize) opts) > 0
                                    then CM.replicateM
                                         ((read zs) * (read (extractBatchSize (DL.head (DL.filter (isBatchSize) opts)))))
                                         (SR.randomRIO (read xs,read ys))
                                    else CM.replicateM
                                         (read zs)
                                         (SR.randomRIO (read xs,read ys))

{----------------------------}

{-Random nucleotides function.-}

--randomNucleotides -> This function will 
--create a list of random integers from 0 
--to 2.
randomNucleotides :: String -> [Flag] -> IO [Int]
randomNucleotides []    []    = return []
randomNucleotides []    (_:_) = return []
randomNucleotides xs    opts  = --Check if user passed BatchSize flag.
                                if DL.length (DL.filter (isBatchSize) opts) > 0
                                    then CM.replicateM 
                                         ((read xs) * (read (extractBatchSize (DL.head (DL.filter (isBatchSize) opts))))) 
                                         (SR.randomRIO (0,2))
                                    else CM.replicateM 
                                         (read xs) 
                                         (SR.randomRIO (0,2))

{------------------------------}


{-Random snv generator function.-}

--randomSnvGenerator -> This function will
--generate random snvs given user input.
randomSnvGenerator :: String -> [Sequence] -> [Int] -> [Int] -> [[String]]
randomSnvGenerator [] [] [] [] = []
randomSnvGenerator as es fs gs = DL.map (\(a,b,c,(d,e)) -> [a,b,c,d,e]) (randomSnvGeneratorSmall as es fs gs) 
--randomSnvGeneratorSmall -> This function will
--generate random snvs given user input.
randomSnvGeneratorSmall :: String -> [Sequence] -> [Int] -> [Int] -> [(String,String,String,(String,String))]
randomSnvGeneratorSmall [] [] [] [] = []
randomSnvGeneratorSmall as es fs gs = randomSnvs as es fs gs
    where
        --Local definitions.--
        --randomSnvs -> This function will
        --generate random positions within
        --specified start and stop.
        randomSnvs :: String -> [Sequence] -> [Int] -> [Int] -> [(String,String,String,(String,String))]
        randomSnvs as es fs gs = randomSnvsSmall (grabFastaSequence (read as) es) fs gs as
        --randomSnvsSmall -> This function will
        --generate random positions within
        --specified start and stop. 
        randomSnvsSmall :: DBC.ByteString -> [Int] -> [Int] -> String -> [(String,String,String,(String,String))]
        randomSnvsSmall _ []    _  _  = []
        randomSnvsSmall _ (_:_) [] _  = []
        randomSnvsSmall xs (y:ys) (z:zs) as = [(as,show y,show y,mapTuple (\x -> [x]) (((DL.concatMap (\s -> DL.filter (\(r,_) -> r == s) nucleotidemapping) [DBC.index xs (y+1)])) DL.!! z))] ++ (randomSnvsSmall xs ys zs as)
        --grabFastaSequence -> This function will
        --grab the correct fasta sequence
        --using chromosome information
        --in the region file.
        grabFastaSequence :: Int -> [Sequence] -> DBC.ByteString
        grabFastaSequence x ys = smallGrabFastaSequence x ys [0..(DL.length ys) - 1]
        --smallGrabFastaSequence -> This function will
        --grab the correct fasta sequence
        --using chromosome information
        --in the region file.
        smallGrabFastaSequence :: Int -> [Sequence] -> [Int] -> DBC.ByteString
        smallGrabFastaSequence _ _ [] = DBC.empty
        smallGrabFastaSequence x ys (z:zs) = if ((bslToStr (extractunSL (extractSeqLabel (ys !! z)))) == ("chr" ++ show x))
                                                 then strToBSC8 (bslToStr (extractunSD (extractSeqData (ys !! z))))
                                                 else smallGrabFastaSequence x ys zs
        --nucleotidemapping -> List containing possible mappings for nucleotides mutations.
        --(NUCLEOTIDE,NUCLEOTIDE_TRANSITION/NUCLEOTIDE_TRANSVERSION)
        nucleotidemapping = [('A','T'),('A','G'),('A','C'),('T','A'),('T','G'),('T','C'),('G','C'),('G','A'),('G','T'),('C','G'),('C','A'),('C','T')]
        ----------------------

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
                           (DL.map (show) 
                           (DL.concat 
                           (DL.map 
                           (\x -> DL.take ((read ys) `div` ((read (extractBatchSize (DL.head (DL.filter (isBatchSize) opts))) :: Int))) 
                           (DL.repeat x)) 
                           [1..(read (extractBatchSize (DL.head (DL.filter (isBatchSize) opts))) :: Int)])))
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
    let tabsandnewlinesadded = DL.map (mapNotLast (++ "\t")) xs
    --Write the output to the user-specified filename.
    SIO.writeFile (outfilestring) $
                  (TPB.render $
                  (TPB.hsep 0 TPB.left . DL.map (TPB.vcat TPB.left) . DL.map (DL.map (TPB.text)))
                  (DL.transpose tabsandnewlinesadded)) 

{---------------------}


{-FRR Specific Function.-}

--processArgsAndFiles -> This function will
--walk through each of the command-line
--arguments and files provided by the user.
processArgsAndFiles :: ([Flag],[String]) -> IO ()
processArgsAndFiles ([],[]) = return () 
processArgsAndFiles (options,files) = do
    --Process the Chromosome argument.
    let readchromosomearg = files DL.!! 0 
    ----------------------------------
    --Process the Start argument.
    let readstartarg = files DL.!! 1
    -----------------------------
    --Process the Stop argument.
    let readstoparg = files DL.!! 2
    ----------------------------
    --Process the Number of batches.
    let readnumbatchesarg = files DL.!! 3
    ---------------------------------------------
    --Process the fasta file.
    readfastafile <- BSF.readFasta (files DL.!! 4)   
    -------------------------
    --Get random positions.
    randompositions <- randomPositions readstartarg readstoparg readnumbatchesarg options
    --Get random nucleotides.
    randomnucleotides <- randomNucleotides readnumbatchesarg options
    --Run randomsnvsgenerator.
    let randomsnvs = randomSnvGenerator readchromosomearg 
                                        readfastafile
                                        randompositions
                                        randomnucleotides
    --Add batch column to randomsnvs.
    let batchadded = batchAdder randomsnvs readnumbatchesarg options 
    --Add header to randomsnvs.
    let finalpreprint = [["Batch","Chromosome","Start","Stop","Ref","Alt"]] ++ batchadded
    --Print the file to stdout (cat) or to a file.
    if DL.length (DL.filter (isOutputFile) options) > 0 
        then printFile options finalpreprint
        else catFile finalpreprint

{-------------------------}


{-Main function.-}

main :: IO ()
main = do
    --Get command line arguments.
    (args,files) <- SE.getArgs >>= compilerOpts
    --Run args and files through processArgsandFiles.
    processArgsAndFiles (args,files)

{----------------}
