#ifndef ROOTGEN_H
#define ROOTGEN_H
#define debug
#include <string>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <utype.h>
#include <utchain.h>

class rootGen
{
public:
    rootGen();
    ~rootGen();
    int stirTemplateGen();
    int createEmptyMichelogram();
    void clearMichelogram();

    void createROOTMichelogram();
    void saveMichelogram();
    int loadROOTfiles();
    void processMacData();
    int showFirstaa();
#ifdef debug
    void debugWork();
#endif

    std::string rootFileName="file.root";
    std::string output_template_name="tst",originating_system="PET",data_size="float",
    byte_order="LITTLEENDIAN",transaxial_sampling;
    //datasize:"short integer" "float" "integer"

    int16_t minimum_ring_difference=0,maximum_ring_difference=99;
    //mac param
    bool has_block,has_module,has_rsector;
    float inner_ring_param1=2.0;//=g_wild->crystal_array[6]
    float inner_ring_param2,inner_ring_param3,inner_ring_param4,inner_ring_param5;//crystal_array[3],block_array[3],module_array.at(3),rsector_array.at(1)
    float average_depth_of_interaction=0,view_offset_degrees=0;
    int16_t arc_corrected_bins=0;
    float image_scaling_factor=1;
    int16_t data_offset=0,time_frames=1;
    int16_t number_of_rings=100;

    int16_t number_of_detector_per_ring=624;
    int16_t number_of_bins=312;
    float inner_ring_diameter=83.5;
    float distance_between_rings=0.2f;

    int16_t blocks_per_bucket_in_tras=1;
    int16_t blocks_per_bucket_in_axial=3;
    int16_t crystals_per_block_axial=13;
    int16_t crystals_per_block_trans=13;
    int16_t detector_layers=1;
    int16_t crystals_per_singles_axial=13;
    int16_t crystals_per_singles_trans=13;
    int16_t suggested_offset;

    float effective_cental_bin_size=0.267628f;

    int16_t matrix_size_4, matrix_size_3, matrix_size_2, matrix_size_1;
    //GATE simulation usefulls
    //useless

    // various types od michelograms
    float **** float_michelogram;
    qint32 **** int_michelogram;
    qint16 **** short_michelogram;

    qint16 selected_ring_difference=100;
    //output settings
    bool save_as_projections=false;

    //params for root event process
    bool load_test=false;
    qint16 detectors_per_ring=960;
    qint16 tang_bins=480;
    std::string scanner_name="PET";
    qint16 crystals_per_module_z=1,  crystals_per_module_xy=25;
    qint16 modules_z=1, modules_xy=8;
    qint16 blocks_z=1, blocks_xy=16;
    qint16 crystals_z=1, crystals_xy=16;
    qint16 offset=0;
    qint16 gate_offset=0;

    float inner_diameter;
    std::string mac_filename="in.mac";

    //Statistics

    unsigned long int trues, scattered, random, prompts, recorded_counts;
    // ROOT tree declarations
    //Event Data

    TChain *Coincidences;

    double nentries=1000;

    Int_t           compton1=1, compton2=1;
    Int_t           runID=1, eventID1=1, eventID2=1;
    Int_t           crystalID1=1, crystalID2=10;
    Int_t           submoduleID1=1, submoduleID2=5, moduleID1=1, moduleID2=4, rsectorID1=1, rsectorID2=2;
    Int_t           rotAngleID1=1, rotAngleID2=10;
    Float_t          energy_1=0.511f, energy_2=0.511f;
    Float_t          globalPosX1=1, globalPosX2=2;
    Int_t           comptonPhantom1=0, comptonPhantom2=0;
    Int_t           sourceID1=1, sourceID2=2;
    Float_t         sourcePosX1=0.1f, sourcePosY1=0.1f, sourcePosZ1=0.1f, sourcePosX2=0.2f, sourcePosY2=0.2f, sourcePosZ2=0.2f;





    //param used in michelogram generation

    int low_energy_u=350,hi_energy_u=650;
    bool include_randoms=false;
    bool include_scattered=false;

















};


#endif // ROOTGEN_H
