from configparser import ConfigParser
import logging
from pathlib import Path
import ROOT
import numpy

import math
from math import floor

log=logging.Logger("Main")
Coincidences = ROOT.TChain("Coincidences")

class iConfigParser(ConfigParser):
    def __init__(self):
        super(iConfigParser,self).__init__(allow_no_value=True,inline_comment_prefixes=';')
        self.optionxform = str

def CreateDefaultConfig()->iConfigParser:
    config=iConfigParser()
    config.add_section("FILE")
    config['FILE']['InputRootFileName']='file.root'
    config['FILE']['OutputNameString']='test'
    config['FILE']['OutputHeaderSuffix']='.hs'
    config['FILE']['OutputDataSuffix']='.s'    
    config['FILE']['OriginatingSystem']='PET'
    config['FILE']['DataType']='float ;datasize: short integer/float/integer'  #datasize:"short integer" "float" "integer"
    config['FILE']['ByteOrder']='LITTLEENDIAN'

    config.add_section("GEOMETRY")
    config['GEOMETRY']
    config['GEOMETRY']['MinimumRingDifference']='0'
    config['GEOMETRY']['MaximumRingDifference']='1'
    config['GEOMETRY']['NumberOfRings']='1'
    config['GEOMETRY']['NumberOfDetectorPerRing']='960'
    config['GEOMETRY']['NumberOfBins']='480'
    config['GEOMETRY']['TangentialBins']='480'
    config['GEOMETRY']['InnerRingDiameter']='83.5 ;in cm'
    config['GEOMETRY']['DistanceBetweenRings']='0.2 ;in cm'
    config['GEOMETRY']['BlocksPerBucketInTrans']='1'
    config['GEOMETRY']['BlocksPerBucketInAxial']='5'   
    config['GEOMETRY']['HasBlock']='0'     
    config['GEOMETRY']['CrystalsPerBlockInTrans']='20'
    config['GEOMETRY']['CrystalsPerBlockInAxial']='20'  
    config['GEOMETRY']['HasModule']='0'
    config['GEOMETRY']['DetectorLayers']='1'
    config['GEOMETRY']['HasRsector']='1'
    config['GEOMETRY']['CrystalsPerSinglesInTrans']='20'
    config['GEOMETRY']['CrystalsPerSinglesInAxial']='20' 
    config['GEOMETRY']['InnerRingParam01']='2.0'
    config['GEOMETRY']['InnerRingParam02']='1.0'
    config['GEOMETRY']['InnerRingParam03']='0.0'
    config['GEOMETRY']['InnerRingParam04']='0.0'
    config['GEOMETRY']['InnerRingParam05']='0.0'
    config['GEOMETRY']['AverageDepthOfInteraction']='0.0'
    config['GEOMETRY']['ViewOffsetDegrees']='0.0'
    config['GEOMETRY']['ImageScalingFactor']='1'
    config['GEOMETRY']['DataOffset']='0'
    config['GEOMETRY']['TimeFrames']='1'
    config.add_section("COINCIDENCE")
    config['COINCIDENCE']['LowEnergy']="350 ;in KeV"
    config['COINCIDENCE']['HighEnergy']="650 ;in KeV"
    config['COINCIDENCE']['IncludeRandoms']="False"
    config['COINCIDENCE']['IncludeScattered']="False"
    config['COINCIDENCE']["LimitRingDifference"]="False"
    config['COINCIDENCE']["SelectedRingDifference"]="1"
    return config


def StirTemplateGen(config:iConfigParser):
        

    OutputFilename=config['FILE']['OutputNameString']+ config['FILE']['OutputHeaderSuffix']

    DataSize=0
    if config['FILE']['DataType']=='short integer':
        DataSize=2
    elif config['FILE']['DataType']== 'float' or config['FILE']['DataType']=='integer':
        DataSize=4
    else:
        DataSize=4
        log.warning("[FILE] DataType not found or illegal, using default DataSize=4")

    NumberOfRings=config['GEOMETRY'].getint('NumberOfRings')
    MatrixSize_4=NumberOfRings*2-1
    NumberOfDetectorPerRing=config['GEOMETRY'].getint('NumberOfDetectorPerRing')
    MatrixSize_3=floor(NumberOfDetectorPerRing/2)
    AxialCoordinate="{"
    if (NumberOfRings %2) !=0:
        for i in range(NumberOfRings):
            AxialCoordinate+=str(i+1)
            AxialCoordinate+=","
        for i in range(NumberOfRings,0,-1):
            AxialCoordinate+=str(i)
            if i>1:
                AxialCoordinate+=","
    else:
        for i in range(NumberOfRings):
            AxialCoordinate+=str(i+1)
            AxialCoordinate+=","
        for i in range(NumberOfRings-1,0,-1):
            AxialCoordinate+=str(i)
            if i>1:
                AxialCoordinate+=","
    AxialCoordinate+=" }"
    MatrixSize_2=AxialCoordinate
    MatrixSize_1=floor(NumberOfDetectorPerRing/2)
    MinRing=-NumberOfRings
    MaxRing=+NumberOfRings
    MinimumRingDifference=config['GEOMETRY'].getint('MinimumRingDifference')
    MaximumRingDifference=config['GEOMETRY'].getint('MaximumRingDifference')
    MinimumRingArray="{"
    if (NumberOfRings %2) !=0:
        for i in range(NumberOfRings-MinimumRingDifference):
            MinimumRingArray+=str(MinRing+(i+1))
            MinimumRingArray+=","
        for i in range(MinimumRingDifference+1,NumberOfRings,+1):
            MinimumRingArray+=str(i)
            if i<(NumberOfRings-1):
                MinimumRingArray+=","
    else:
        for i in range(NumberOfRings-MinimumRingDifference):
            MinimumRingArray+=str(MinRing+(i+1))
            MinimumRingArray+=","
        for i in range(MinimumRingDifference+1,NumberOfRings,+1):
            MinimumRingArray+=str(i)
            if i<(NumberOfRings-1):
                MinimumRingArray+=","
    MinimumRingArray+="}"


    MaximumRingArray="{"
    if (NumberOfRings %2) !=0:
        for i in range(-MaximumRingDifference,0,1):
            MaximumRingArray+=str(i)
            MaximumRingArray+=","
        for i in range(0,MaximumRingDifference+1,+1):
            MaximumRingArray+=str(i)
            if i<(MaximumRingDifference):
                MaximumRingArray+=","
    else:
        for i in range(-MaximumRingDifference,0,1):
            MaximumRingArray+=str(i)
            MaximumRingArray+=","
        for i in range(0,MaximumRingDifference+1,+1):
            MaximumRingArray+=str(i)
            if i<(MaximumRingDifference):
                MaximumRingArray+=","
    MaximumRingArray+="}"


    InnerRingParam01=config['GEOMETRY'].getfloat('InnerRingParam01')
    InnerRingParam02=config['GEOMETRY'].getfloat('InnerRingParam02')
    InnerRingParam03=config['GEOMETRY'].getfloat('InnerRingParam03')
    InnerRingParam04=config['GEOMETRY'].getfloat('InnerRingParam04') 
    InnerRingParam05=config['GEOMETRY'].getfloat('InnerRingParam05')   
    HasBlock=config['GEOMETRY'].getboolean('HasBlock')     
    HasModule=config['GEOMETRY'].getboolean('HasModule')      
    HasRsector=config['GEOMETRY'].getboolean('HasRsector')    
    tInnerRingDiameter=-InnerRingParam01/2.0+InnerRingParam02
    if HasBlock:
        tInnerRingDiameter+=InnerRingParam03
    if HasModule:
        tInnerRingDiameter+=InnerRingParam04
    if HasRsector:
        tInnerRingDiameter+=InnerRingParam05
    tInnerRingDiameter*=2.0
    EffectiveCentralBinSize=tInnerRingDiameter/NumberOfRings


    HeaderTemplate="""    !INTERFILE :=
    name of data file := {OutputFilename}
    originating system := {OriginatingSystem}
    !GENERAL DATA :=
    !GENERAL IMAGE DATA :=
    !type of data := PET 
    imagedata byte order := {ByteOrder}
    !PET STUDY (General) := 
    !PET data type := Emission
    applied corrections := {{None}}
    !number format := {DataTypeName}
    !number of bytes per pixel := {DataSize}
    number of dimensions := 4
    matrix axis label [4] := segment
    !matrix size [4] := {MatrixSize_4}
    matrix axis label [3] := view
    !matrix size [3] := {MatrixSize_3}
    matrix axis label [2] := axial coordinate
    !matrix size [2] := {MatrixSize_2}
    matrix axis label [1] := tangential coordinate
    !matrix size [1] := {MatrixSize_1}
    minimum ring difference per segment := {MinimumRingDifference}
    maximum ring difference per segment := {MaximumRingDifference}
    Scanner parameters:= 
    Scanner type := {ScannerType}
    Number of rings                          := {NumberOfRings}
    Number of detectors per ring             := {NumberOfDetectorPerRing}
    Inner ring diameter (cm)                 := {InnerRingDiameter}
    Average depth of interaction (cm)        := {AverageDepthofInteraction}
    Distance between rings (cm)              := {DistanceBetweenRings}
    Default bin size (cm)                    := {DefaultBinSize}
    View offset (degrees)                    := {ViewOffset}
    Maximum number of non-arc-corrected bins := {MaxNumofNACBins}
    Default number of arc-corrected bins     := {DefaultNumofACBins}
    Number of blocks per bucket in transaxial direction         := {BlocksPerBucketInTrans}
    Number of blocks per bucket in axial direction              := {BlocksPerBucketInAxial}
    Number of crystals per block in transaxial direction        := {CrystalsPerBlockInTrans}
    Number of crystals per block in axial direction             := {CrystalsPerBlockInAxial}
    Number of detector layers                                   := {DetectorLayers}
    Number of crystals per singles unit in transaxial direction := {CrystalsPerSinglesInTrans}
    Number of crystals per singles unit in axial direction      := {CrystalsPerSinglesInAxial}
    end scanner parameters:=
    effective central bin size (cm) := {EffectiveCentralBinSize}
    image scaling factor[1] := {ImageScalingFactor}
    data offset in bytes[1] := {DataOffset}
    number of time frames := {TimeFrames}
    !END OF INTERFILE :=
    """


    InterVarDict={
        "OutputFilename":OutputFilename,
        "OriginatingSystem":config['FILE']["OriginatingSystem"],
        "ByteOrder":config['FILE']["ByteOrder"],
        "DataTypeName":config['FILE']['DataType'],
        "DataSize":DataSize,
        "MatrixSize_4":MatrixSize_4,
        "MatrixSize_3":MatrixSize_3,
        "MatrixSize_2":MatrixSize_2,
        "MatrixSize_1":MatrixSize_1,
        "MinimumRingDifference":config['GEOMETRY'].getint('MinimumRingDifference'),
        "MaximumRingDifference":config['GEOMETRY'].getint('MaximumRingDifference'),
        "ScannerType":config['FILE']["OriginatingSystem"],
        "NumberOfRings":config['GEOMETRY']['NumberOfRings'],
        "NumberOfDetectorPerRing":config['GEOMETRY']['NumberOfDetectorPerRing'],
        "InnerRingDiameter":tInnerRingDiameter,
        "AverageDepthofInteraction":0,
        "DistanceBetweenRings":config['GEOMETRY']['DistanceBetweenRings'],
        "DefaultBinSize":config['GEOMETRY']['DistanceBetweenRings'],
        "ViewOffset":config['GEOMETRY']['ViewOffsetDegrees'],
        "MaxNumofNACBins":config['GEOMETRY']['NumberOfBins'],
        "DefaultNumofACBins":config['GEOMETRY']['TangentialBins'],
        "BlocksPerBucketInTrans":config['GEOMETRY']['BlocksPerBucketInTrans'],
        "BlocksPerBucketInAxial":config['GEOMETRY']['BlocksPerBucketInAxial'],       
        "CrystalsPerBlockInTrans":config['GEOMETRY']['CrystalsPerBlockInTrans'],
        "CrystalsPerBlockInAxial":config['GEOMETRY']['CrystalsPerBlockInAxial'],  
        "DetectorLayers":config['GEOMETRY']['DetectorLayers'],
        "CrystalsPerSinglesInTrans":config['GEOMETRY']['CrystalsPerSinglesInTrans'],
        "CrystalsPerSinglesInAxial":config['GEOMETRY']['CrystalsPerSinglesInAxial'], 
        "EffectiveCentralBinSize":EffectiveCentralBinSize,
        "ImageScalingFactor":config['GEOMETRY'].getfloat('ImageScalingFactor'),
        "DataOffset":config['GEOMETRY'].getfloat('DataOffset'),
        "TimeFrames":config['GEOMETRY'].getint('TimeFrames'),
        "LowEnergy":config['COINCIDENCE']['LowEnergy'],
        "HighEnergy":config['COINCIDENCE']['HighEnergy']
    }
    InterVarValidationCheck(InterVarDict)
    Header=HeaderTemplate.format_map(InterVarDict)
    HeaderFilename=config['FILE']['OutputNameString']+config['FILE']['OutputHeaderSuffix']
    with open(HeaderFilename,"w") as f:
        f.write(Header)
    return 1


def InterVarValidationCheck(InterVar:dict)->bool:
    print(InterVar["DataTypeName"])
    assert (InterVar["DataTypeName"] in ["short integer","float","integer"])==True
    return True




def CreateROOTMichelogram(iMichelogram:numpy.array,config:iConfigParser):

    infileName=config['FILE']['InputRootFileName']
    # Coincidences.Add(infileName)
    # Coincidences.SetBranchStatus("*",0) ## NOT working for PyROOT
    inFile=ROOT.TFile.Open("file.root","READ")
    Coincidences=inFile.Get("Coincidences")
    entries=Coincidences.GetEntries()
    print("entries:",entries)


    IncludeRandoms=config['COINCIDENCE'].getboolean('IncludeRandoms')
    IncludeScattered=config['COINCIDENCE'].getboolean('IncludeScattered')
    LimitRingDifference=config['COINCIDENCE'].getboolean('LimitRingDifference')
    SelectedRingDifference=1
    if LimitRingDifference:
        SelectedRingDifference=config['COINCIDENCE'].getint('SelectedRingDifference')
    else:
        SelectedRingDifference=config['GEOMETRY'].getint('MaximumRingDifference')-config['GEOMETRY'].getint('MinimumRingDifference')+1
    
    LowEnergy = config['COINCIDENCE'].getfloat("LowEnergy")  #in KeV
    LowEnergy /=1000 #in MeV
    HighEnergy = config['COINCIDENCE'].getfloat("HighEnergy") #in KeV
    HighEnergy /=1000 #in MeV
    InnerRingDiameter =config['GEOMETRY'].getfloat("InnerRingDiameter")
    DetectorPerRing=config['GEOMETRY'].getint('NumberOfDetectorPerRing')
    TangentialBins=config['GEOMETRY'].getint('TangentialBins')
    crystals_xy=config["GEOMETRY"].getint("CrystalsPerBlockInTrans")
    crystals_z=config["GEOMETRY"].getint("CrystalsPerBlockInAxial")


    blocks_xy=1 #BlockPerModuleInTrans Block:=SubModule
    blocks_z=1 #BlockPerModuleInAxial Block:=SubModule
    modules_xy=1 #ModulePerRingInTrans
    modules_z=5 #ModulePerRingInAxial
    
    Trues = 0
    Scattered = 0
    Random = 0
    Recorded_counts = 0
    Prompts = 0
    AverageDepthofInteraction = 0
    for entryNum in range(0,Coincidences.GetEntries()):
        Coincidences.GetEntry(entryNum)
        energy_1=Coincidences.energy1
        print(energy_1)
        energy_2=Coincidences.energy2
        if energy_1<LowEnergy or energy_1>HighEnergy or energy_2<LowEnergy or energy_2>HighEnergy:
            continue
        Prompts+=1
        
        eventID1=Coincidences.eventID1
        eventID2=Coincidences.eventID2
        comptonPhantom1=Coincidences.comptonPhantom1
        comptonPhantom2=Coincidences.comptonPhantom2
        if eventID1==eventID2:
            if (comptonPhantom1 == 0 and comptonPhantom2 ==0):
                Trues+=1
            else:
                Scattered+=1
                if not IncludeScattered:
                    continue
        
        else:
            Random+=1
            if not IncludeRandoms:
                continue
        
        crystalID1=Coincidences.crystalID1
        crystalID2=Coincidences.crystalID2
        submoduleID1=Coincidences.submoduleID1
        submoduleID2=Coincidences.submoduleID2
        moduleID1=Coincidences.moduleID1
        moduleID2=Coincidences.moduleID2
        rsectorID1=Coincidences.rsectorID1
        rsectorID2=Coincidences.rsectorID2

        ### IDENTIFY THE RING  
        ring1=crystalID1//crystals_xy
        ring1+=(submoduleID1//blocks_xy)*crystals_z
        ring1+=(moduleID1// modules_xy)*blocks_z*crystals_z

        ring2=crystalID2//crystals_xy
        ring2+=(submoduleID2//blocks_xy)*crystals_z
        ring2+=(moduleID2// modules_xy)*blocks_z*crystals_z

        if abs(ring1-ring2)>SelectedRingDifference:
            continue  
            
        

        ### IDENTIFY THE CRYSTAL
        crystal1 = rsectorID1 *  modules_xy *  blocks_xy *  crystals_xy \
                + (moduleID1% modules_xy) *  blocks_xy *  crystals_xy \
                + (submoduleID1% blocks_xy) *  crystals_xy \
                + (crystalID1% crystals_xy)
        crystal2 = rsectorID2 *  modules_xy *  blocks_xy *  crystals_xy \
                + (moduleID2% modules_xy) *  blocks_xy *  crystals_xy \
                + (submoduleID2% blocks_xy) *  crystals_xy \
                + (crystalID2% crystals_xy)

        ### Rotate the image
        
        GateOffset=floor(DetectorPerRing*0.75) ### Rotate 90 degrees
        
        if GateOffset>0:
            crystal1 = crystal1 +  GateOffset
            crystal2 = crystal2 +  GateOffset
            if (crystal1 >=  DetectorPerRing):
                crystal1 = crystal1 -  DetectorPerRing
            if (crystal2 >=  DetectorPerRing):
                crystal2 = crystal2 -  DetectorPerRing
        
        ##--------------------------------
        ##  Bin the crystal ring pairs into Michelograms
        ##  u - radial sinogram component
        ##  phi - azimuthal sinogram component
        ##  ring pairs are sorted according to c1 < c2 else flip
        ##  where c1 and c2 are crystals at phi(u = S_WIDTH/2)
        ##--------------------------------

        phi= ((crystal1 + crystal2 +  DetectorPerRing//2)% DetectorPerRing)//2
        if (((crystal1 + crystal2) < (3* DetectorPerRing//2)) and ((crystal1 + crystal2) >= ( DetectorPerRing//2))):
            u    =  abs(crystal1 - crystal2) -   DetectorPerRing//2 +  TangentialBins//2
        else:
            u = -abs(crystal1 - crystal2) +   DetectorPerRing//2 +  TangentialBins//2

        if ( u >=  TangentialBins or u < 0 ):
            continue
        if (u%2 == 0):
            zi = ( DetectorPerRing//2 - (crystal1 - crystal2) - 1)//2
            if (zi >=   DetectorPerRing//4):
                zi = zi -  DetectorPerRing//2 + 1
            if (zi <= - DetectorPerRing//4):
                zi = zi +  DetectorPerRing//2 - 1

        else:
            zi = ( DetectorPerRing//2 - (crystal1 - crystal2))//2
            if (zi >=   DetectorPerRing//4):
                zi = zi -  DetectorPerRing//2
            if (zi <= - DetectorPerRing//4):
                zi = zi +  DetectorPerRing//2

        c1 = crystal1 + zi
        c2 = crystal2 - zi
        if (c1 >=  DetectorPerRing):
            c1 = c1 -  DetectorPerRing
        elif (c1 < 0):
            c1 = c1 +  DetectorPerRing
        if (c2 >=  DetectorPerRing):
            c2 = c2 -  DetectorPerRing
        elif (c2 < 0):
            c2 = c2 +  DetectorPerRing
        if (c1 < c2):
            flip = 0
        else:
            flip = 1
        if (flip):
            swap  = ring1
            ring1 = ring2
            ring2 = swap

        iMichelogram[ring2][ring1][phi][u]+=1.0

        Recorded_counts+=1
        #AverageDepthofInteraction+= 0 # Have not Implemented

    #AverageDepthofInteraction /= Recorded_counts
    #AverageDepthofInteraction -= InnerRingDiameter/2.0
    print("Trues=%d" % Trues)
    print("Scattered=%d" % Scattered)
    print("Random=%d" % Random)
    print("Recorded_counts=%d" % Recorded_counts)
    print("Prompts=%d" % Prompts)
    
    
    
    return iMichelogram
def GetSelectedRingDifference(config:iConfigParser):
    LimitRingDifference=config['COINCIDENCE'].getboolean('LimitRingDifference')
    SelectedRingDifference=1
    if LimitRingDifference:
        SelectedRingDifference=config['COINCIDENCE'].getint('SelectedRingDifference')
    else:
        SelectedRingDifference=config['GEOMETRY'].getint('MaximumRingDifference')-config['GEOMETRY'].getint('MinimumRingDifference')+1
    return SelectedRingDifference

def CreateEmptyFloatMichelogram(config:iConfigParser):
    NumberOfRings=config['GEOMETRY'].getint('NumberOfRings')
    DetectorPerRing=config['GEOMETRY'].getint('NumberOfDetectorPerRing')
    TangentialBins=config['GEOMETRY'].getint('TangentialBins')    
    m1=NumberOfRings
    m2=NumberOfRings
    m3=floor(DetectorPerRing//2 +1.5)
    m4=TangentialBins
    mm=GetSelectedRingDifference(config)
    assert mm>=m1
    iMichelogram=numpy.zeros((m1,m2,m3,m4),dtype=float)
    return iMichelogram

def main():
    
    ConfPath=Path('conf.ini')
    if not ConfPath.exists():
        config_d=CreateDefaultConfig()
        with open(ConfPath,'w') as f:
            config_d.write(f)
    config=iConfigParser()
    with open(ConfPath,'r') as f:
        config.read_file(f)
    
    StirTemplateGen(config)
    iMichelogram=CreateEmptyFloatMichelogram(config)
    CreateROOTMichelogram(iMichelogram,config)

    SaveROOTMichelogram(iMichelogram,config)

    
def SaveROOTMichelogram(nparray:numpy.ndarray,config):
    OutputNameString=config['FILE']['OutputNameString']
    OutputDataSuffix=config['FILE']['OutputDataSuffix']
    dst=Path(OutputNameString+OutputDataSuffix)
    nparray.tofile(dst)
    return 1
    
if __name__ =='__main__':
    main()