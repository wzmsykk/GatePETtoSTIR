#include "rootgen.h"

rootGen::rootGen()
{
    Coincidences = new TChain("Coincidences");
}
rootGen::~rootGen()
{

    delete Coincidences;
}
int rootGen::stirTemplateGen()
{
    //std::string *output_template_name, std::string *_originating_system, std::string *_data_size,
    //std::string *_byte_order, bool *_transaxial_sampling,
    //int16_t *_minimum_ring_difference, int16_t *_maximum_ring_difference,
    //qreal *_average_depth_of_interaction, qreal *_view_offset_degrees,
    //int16_t *_arc_corrected_bins, qreal *_image_scaling_factor,
    //int16_t *_data_offset, int16_t *_time_frames)

        std::string output_data_name = output_template_name;
        std::string output_header_name = output_template_name+".hs";
        output_data_name = output_data_name+".s";
        std::string matrix_size, view_size, axial_coordinate, minimum_ring_array, maximum_ring_array;

        int16_t min_ring = - number_of_rings;
        int16_t max_ring = number_of_rings;

        std::string temp;
        //int16_t index_to_write = 0;
        std::string temp_2;



        //param processing
        if (data_size.find("short integer")!=std::string::npos)
            temp_2 = '2';
        else if (data_size.find("float")!=std::string::npos || data_size.find("integer")!=std::string::npos)
            temp_2 = '4';
        else temp_2 = '4';
        matrix_size_4 = (number_of_rings*2)-1;
        std::string size_4_out = std::to_string(matrix_size_4);
        matrix_size_3 = number_of_detector_per_ring/2;
        std::string size_3_out = std::to_string(matrix_size_3);
        axial_coordinate = "{ ";
//                    if (*_transaxial_sampling)
        {
            if (number_of_rings % 2 !=0)
            {
                for (int16_t j = 0; j < number_of_rings;j++)
                {
                    axial_coordinate.append(std::to_string(j+1));
                    axial_coordinate.append(",");
                }

                for (int16_t j = number_of_rings-1; j >0 ; j--)
                {
                    axial_coordinate.append(std::to_string(j));
                    if (j>1)
                        axial_coordinate.append(",");
                }
            }
            else
            {
                for (int16_t j = 0; j < number_of_rings;j++)
                {
                    axial_coordinate.append(std::to_string(j+1));
                    axial_coordinate.append(",");
                }

                for (int16_t j = number_of_rings-1; j >0 ; j--)
                {
                    axial_coordinate.append(std::to_string(j));
                    if (j>1)
                        axial_coordinate.append(",");
                }
            }
        }
        axial_coordinate.append(" }");
        std::string size_2_out=axial_coordinate;
        matrix_size = std::to_string(number_of_rings);
        std::string size_1_out = std::to_string(number_of_detector_per_ring/2);
        minimum_ring_array = "{ ";
        if (number_of_rings%2 !=0)
        {
            for (int16_t j = 0; j < (number_of_rings - minimum_ring_difference);j++)
            {
                minimum_ring_array.append(std::to_string(min_ring+(j+1)));
                minimum_ring_array.append(",");
            }

            for (int16_t j = minimum_ring_difference+1; j < number_of_rings ; j++)
            {
                minimum_ring_array.append(std::to_string(j));
                if (j<number_of_rings - 1)
                    minimum_ring_array.append(",");
            }
        }
        else
        {
            for (int16_t j = 0; j < (number_of_rings - minimum_ring_difference);j++)
            {
                minimum_ring_array.append(std::to_string(min_ring+(j+1)));
                minimum_ring_array.append(",");
            }

            for (int16_t j = minimum_ring_difference+1; j < number_of_rings ; j++)
            {
                minimum_ring_array.append(std::to_string(j));
                if (j<number_of_rings -1)
                    minimum_ring_array.append(",");
            }
        }
        minimum_ring_array.append(" }");

        matrix_size = std::to_string(number_of_rings);
        maximum_ring_array = "{ ";

        if (number_of_rings%2 !=0)
        {
            for (int16_t j = -maximum_ring_difference; j < 0 ; j++)
            {
                maximum_ring_array.append(std::to_string(j));
                //if (j<number_of_rings - 1)
                maximum_ring_array.append(",");
            }


            for (int16_t j = 0; j <= maximum_ring_difference;j++)
            {
                maximum_ring_array.append(std::to_string(j));
                if (j<maximum_ring_difference)
                    maximum_ring_array.append(",");
            }
        }
        else
        {
            for (int16_t j = -maximum_ring_difference; j < 0 ; j++)
            {
                maximum_ring_array.append(std::to_string(j));
                //if (j<number_of_rings - 1)
                maximum_ring_array.append(",");
            }


            for (int16_t j = 0; j <= maximum_ring_difference;j++)
            {
                maximum_ring_array.append(std::to_string(j));
                if (j<maximum_ring_difference)
                    maximum_ring_array.append(",");
            }
        }

        maximum_ring_array.append(" }");

        matrix_size = std::to_string(number_of_rings);

        arc_corrected_bins=tang_bins;
        inner_ring_diameter = -inner_ring_param1/2.0;
            inner_ring_diameter = inner_ring_diameter + inner_ring_param2;
            if (has_block)
                inner_ring_diameter = inner_ring_diameter + inner_ring_param3;
            if (has_module)
                inner_ring_diameter = inner_ring_diameter + inner_ring_param4;
            if (has_rsector)
                inner_ring_diameter = inner_ring_diameter + inner_ring_param5;
            inner_ring_diameter = inner_ring_diameter * 2.0;
        blocks_per_bucket_in_tras=modules_xy;
        blocks_per_bucket_in_axial=modules_z;
        crystals_per_block_trans=crystals_xy;
        crystals_per_block_axial=crystals_z;
        crystals_per_singles_trans=crystals_xy;
        crystals_per_singles_axial=crystals_z;
        effective_cental_bin_size=inner_ring_diameter / number_of_bins;
        image_scaling_factor=1;
        time_frames=1;
        data_offset=0;
//writefiles
        std::ofstream out(output_header_name);
        out <<"!INTERFILE :="<<std::endl;
        out << "name of data file := "<< output_data_name <<std::endl;
        out << "originating system := "<<originating_system<<std::endl;
        out << "!GENERAL DATA :="<<std::endl;
        out << "!GENERAL IMAGE DATA :="<<std::endl;
        out << "!type of data := PET "<<std::endl;
        out << "imagedata byte order := "<<byte_order<<std::endl;
        out << "!PET STUDY (General) := "<<std::endl;
        out << "!PET data type := Emission"<<std::endl;
        out << "applied corrections := {None}"<<std::endl;
        out << "!number format := "<<data_size<<std::endl;
        out << "!number of bytes per pixel := "<<temp_2<<std::endl;
        out << "number of dimensions := 4"<<std::endl;
        out << "matrix axis label [4] := segment"<<std::endl;
        out << "!matrix size [4] := "<<size_4_out<<std::endl;
        out << "matrix axis label [3] := view"<<std::endl;
        out << "!matrix size [3] := "<<size_3_out<<std::endl;
        out << "matrix axis label [2] := axial coordinate"<<std::endl;
        out << "!matrix size [2] := "<<size_2_out<<std::endl;
        out << "matrix axis label [1] := tangential coordinate"<<std::endl;
        out << "!matrix size [1] := "<<size_1_out<<std::endl;
        out << "minimum ring difference per segment := "<<minimum_ring_array<<std::endl;
        out << "maximum ring difference per segment := "<<maximum_ring_array<<std::endl;
        out << "Scanner parameters:= "<<std::endl;
        out << "Scanner type := "<<originating_system<<std::endl;
        out << "Number of rings                          := "<<number_of_rings<<std::endl;
        out << "Number of detectors per ring             := "<<number_of_detector_per_ring<<std::endl;
        out << "Inner ring diameter (cm)                 := "<<inner_ring_diameter<<std::endl;
        out << "Average depth of interaction (cm)        := "<<average_depth_of_interaction<<std::endl; // Θα το πέρνει από το gate είτε από το πρότυπο??
        out << "Distance between rings (cm)              := "<<distance_between_rings<<std::endl;
        out << "Default bin size (cm)                    := "<<distance_between_rings<<std::endl;
        out << "View offset (degrees)                    := "<<view_offset_degrees<<std::endl;
        out << "Maximum number of non-arc-corrected bins := "<<number_of_bins<<std::endl;
        out << "Default number of arc-corrected bins     := "<<arc_corrected_bins<<std::endl;
        out << "Number of blocks per bucket in transaxial direction         := "<<blocks_per_bucket_in_tras<<std::endl;
        out << "Number of blocks per bucket in axial direction              := "<<blocks_per_bucket_in_axial<<std::endl;
        out << "Number of crystals per block in axial direction             := "<<crystals_per_block_axial<<std::endl;
        out << "Number of crystals per block in transaxial direction        := "<<crystals_per_block_trans<<std::endl;
        out << "Number of detector layers                                   := "<<detector_layers<<std::endl;
        out << "Number of crystals per singles unit in axial direction      := "<<crystals_per_singles_axial<<std::endl;
        out << "Number of crystals per singles unit in transaxial direction := "<<crystals_per_singles_trans<<std::endl;
        out << "end scanner parameters:="<<std::endl;
        out << "effective central bin size (cm) := "<<effective_cental_bin_size<<std::endl;
        out << "image scaling factor[1] := "<<image_scaling_factor<<std::endl;
        out << "data offset in bytes[1] := "<<data_offset<<std::endl;
        out << "number of time frames := "<<time_frames<<std::endl;
        out << "!END OF INTERFILE :="<<std::endl;

        out.close();
        return 1;
}

int rootGen::createEmptyMichelogram()
{
    if (data_size=="float")
    {

        int16_t m1,m2,m3,m4;
        m1=m2=selected_ring_difference;
        m3=static_cast<int16_t>((detectors_per_ring/2));
        m4=static_cast<int16_t>(tang_bins);
        int16_t mm=static_cast<int16_t>(selected_ring_difference);
        //int16_t mm=m1;
        std::cout<<"Martix params"<<std::endl<<"d1:"<<mm<<" d2:"<<mm<<" d3:"<<m3<<" d4:"<<m4<<std::endl;
        if(mm>m1) {std::cout <<"ERROR"; return 0;}
        std::cout<<"malloc process"<<std::endl;
        float process;
        float_michelogram = static_cast<float****> ( malloc (static_cast<unsigned long>(mm) * sizeof(float***)));

        for (int16_t counter_rings_1 = 0; counter_rings_1< mm; counter_rings_1++)
        {
            process=static_cast<float>(counter_rings_1)/mm;
             std::cout<<'\r'<<"malloc finished %"<<process*100+1;
             std::cout.flush();
            float_michelogram[counter_rings_1] = static_cast<float***>( malloc ( static_cast<unsigned long>(mm) * sizeof(float**)));

            for (int16_t counter_rings_2 = 0; counter_rings_2< mm; counter_rings_2++)
            {

                float_michelogram[counter_rings_1][counter_rings_2] = static_cast<float**>( malloc (static_cast<unsigned long>(m3) * sizeof(float*)));

                for (int16_t counter_dets = 0; counter_dets < m3; counter_dets++)
                {
                    //std::cout<<counter_dets<<std::endl;
                    float_michelogram[counter_rings_1][counter_rings_2][counter_dets] =static_cast<float*>(malloc (static_cast<unsigned long>(m4) * sizeof(float)));
                }
            }
        }
        std::cout<<std::endl;
        std::cout<<"Martix gen OK"<<std::endl;

    }

    else if (data_size=="integer")
    {
        int_michelogram = static_cast<int32_t****>(  malloc ( static_cast<float>(number_of_rings) * sizeof(int32_t***)));

        for (int16_t counter_rings_1 = 0; counter_rings_1< selected_ring_difference; counter_rings_1++)
        {
            int_michelogram[counter_rings_1] = static_cast<int32_t***>( malloc ( static_cast<float>(number_of_rings) * sizeof(int32_t**)));

            for (int16_t counter_rings_2 = 0; counter_rings_2< selected_ring_difference; counter_rings_2++)
            {
                int_michelogram[counter_rings_1][counter_rings_2] = static_cast<int32_t**>( malloc (static_cast<float>(static_cast<int16_t>((detectors_per_ring/2) + 1.5)) * sizeof(int32_t*)));

                for (int16_t counter_dets = 0; counter_dets < static_cast<int16_t>((detectors_per_ring/2) + 1.5); counter_dets++)
                {
                    int_michelogram[counter_rings_1][counter_rings_2][counter_dets] =static_cast<int32_t*>( malloc (static_cast<float>(tang_bins) * sizeof(int32_t)));
                }
            }
        }


    }

    else if (data_size=="short integer")
    {
        short_michelogram = static_cast<int16_t****>(  malloc ( static_cast<float>(number_of_rings) * sizeof(int16_t***)));

        for (int16_t counter_rings_1 = 0; counter_rings_1< selected_ring_difference; counter_rings_1++)
        {
            short_michelogram[counter_rings_1] = static_cast<int16_t***>( malloc ( static_cast<float>(number_of_rings) * sizeof(int16_t**)));

            for (int16_t counter_rings_2 = 0; counter_rings_2< selected_ring_difference; counter_rings_2++)
            {
                short_michelogram[counter_rings_1][counter_rings_2] = static_cast<int16_t**>( malloc (static_cast<float>(static_cast<int16_t>((detectors_per_ring/2) + 1.5)) * sizeof(int16_t*)));

                for (int16_t counter_dets = 0; counter_dets < static_cast<int16_t>((detectors_per_ring/2) + 1.5); counter_dets++)
                {
                    short_michelogram[counter_rings_1][counter_rings_2][counter_dets] =static_cast<int16_t*>( malloc (static_cast<float>(tang_bins) * sizeof(int16_t)));
                }
            }
        }

    }

    return 1;
}



void rootGen::clearMichelogram()
{
    std::cout<<"clearing"<<std::endl;
    float process;
    if (data_size=="float")
    {
        for (qint16 i = 0; i< selected_ring_difference;i++)
        {
            process=static_cast<float>(i)/selected_ring_difference;
            std::cout<<'\r'<<"cleared %"<<process*100+1;
            std::cout.flush();
            for(qint16 j=0; j< selected_ring_difference;j++)
            {
                for(qint16 w=0; w< static_cast<qint16>((detectors_per_ring/2));w++)
                {
                    for (qint16 l=0; l< tang_bins; l++)
                    {
                        float_michelogram[i][j][w][l]=0.0;
                    }
                }
            }
        }
        std::cout<<std::endl;
    }
    else if (data_size=="integer")
    {
        for (qint16 i = 0; i< selected_ring_difference;i++)
        {
            for(qint16 j=0; j< selected_ring_difference;j++)
            {
                for(qint16 w=0; w< static_cast<qint16>((detectors_per_ring/2) + 1.5);w++)
                {
                    for (qint16 l=0; l<tang_bins; l++)
                    {
                        int_michelogram[i][j][w][l]=0;
                    }
                }
            }
        }

    }
    else if (data_size=="short integer")
    {

        for (qint16 i = 0; i< selected_ring_difference;i++)
        {
            for(qint16 j=0; j< selected_ring_difference;j++)
            {
                for(qint16 w=0; w< static_cast<qint16>((detectors_per_ring/2) + 1.5);w++)
                {
                    for (qint16 l=0; l<tang_bins; l++)
                    {
                        short_michelogram[i][j][w][l]=0;
                    }
                }
            }
        }
    }
}




void rootGen::createROOTMichelogram()
{

    //needed params detectors_per_ring crystals_xy
#ifdef debug
  debugWork();
#endif
    Int_t ring1, ring2;
    Int_t crystal1, crystal2;
    std::cout << "MichelogramCreating" << std::endl;
    trues = 0;
    scattered = 0;
    random = 0;
    recorded_counts = 0;
    prompts = 0;
    average_depth_of_interaction = 0;

    float low_energy = static_cast<float>(low_energy_u)/1000;
    float hi_energy =  static_cast<float>(hi_energy_u)/1000;

    //bool include_randoms = ui->checkBox->isChecked();
    //bool include_scattered = ui->checkBox_2->isChecked();
    //defined in rootgen class

    Int_t phi , u;
    Int_t flip, swap, zi, c1, c2;

    // gate_offset = ui->lineEdit_15->text().toShort(&ok, 10);
    //defined in rootgen class


    //prompts_from_each_source.clear();
    //prompts_from_each_source.empty();

    for (unsigned long int i = 0; i < nentries; i++)
    {
        if(i<10){
            std::cout<<"ev1="<<eventID1<<" ev2="<<eventID2<<" comp1="<<comptonPhantom1<<std::endl;
            std::cout<<"cid1="<<crystalID1<<" cid2="<<crystalID2<<" mID1="<<moduleID1<<std::endl;
            std::cout<<"rid1="<<rsectorID1<<" rid2="<<rsectorID2<<" mID2="<<moduleID2<<std::endl;
        }

        //pdialog.setValue(i/nentries);
        //pdialog.show();
        Coincidences->GetEntry(i);

        if (energy_1<low_energy || energy_1>hi_energy || energy_2<low_energy || energy_2>hi_energy)
            continue;

        prompts++;


        if (eventID1 == eventID2)
        {
            if (comptonPhantom1 == 0 && comptonPhantom2 ==0)
            {
                trues++;
                //std::cout << "TrueGot" << std::endl;
            }
            else
            {
                ////Scattered.
                scattered++;
                if (!include_scattered)
                    continue;
            }
        }
        else
        {
            //// Random.
            random++;
            if (!include_randoms)
                continue;
        }


        // Energy Spectra population




        //-----------------------------------
        //  Identify the ring#...
        //-----------------------
        ring1 = static_cast<Int_t>((crystalID1/crystals_xy))
                + static_cast<Int_t>((submoduleID1/ blocks_xy))*crystals_z
                + static_cast<Int_t>((moduleID1/ modules_xy))*blocks_z*crystals_z;
        ring2 = static_cast<Int_t>((crystalID2/crystals_xy))
                + static_cast<Int_t>((submoduleID2/ blocks_xy))* crystals_z
                + static_cast<Int_t>((moduleID2/ modules_xy))* blocks_z* crystals_z;


        if (abs(ring1 - ring2)>selected_ring_difference)
            continue;

        //-----------------------
        //  Identify the crystal#...
        //-----------------------------------
        crystal1 = rsectorID1 *  modules_xy *  blocks_xy *  crystals_xy
                + (moduleID1% modules_xy) *  blocks_xy *  crystals_xy
                + (submoduleID1% blocks_xy) *  crystals_xy
                + (crystalID1% crystals_xy);
        crystal2 = rsectorID2 *  modules_xy *  blocks_xy *  crystals_xy
                + (moduleID2% modules_xy) *  blocks_xy *  crystals_xy
                + (submoduleID2% blocks_xy) *  crystals_xy
                + (crystalID2% crystals_xy);

        //-----------------------------------------------------
        //  Rotate the image correctly#...
        //--------------------------------
        if ( gate_offset > 0)
        {
            crystal1 = crystal1 +  gate_offset;
            crystal2 = crystal2 +  gate_offset;
            if (crystal1 >=  detectors_per_ring)
                crystal1 = crystal1 -  detectors_per_ring;
            if (crystal2 >=  detectors_per_ring)
                crystal2 = crystal2 -  detectors_per_ring;
        }

        //--------------------------------
        //  Bin the crystal ring pairs into Michelograms
        //  u - radial sinogram component
        //  phi - azimuthal sinogram component
        //  ring pairs are sorted according to c1 < c2 else flip
        //  where c1 and c2 are crystals at phi(u = S_WIDTH/2)
        //--------------------------------
        phi = ((crystal1 + crystal2 +  detectors_per_ring/2)% detectors_per_ring)/2;

        if (((crystal1 + crystal2) < (3* detectors_per_ring/2)) && ((crystal1 + crystal2) >= ( detectors_per_ring/2)))
            u    =  abs(crystal1 - crystal2) -   detectors_per_ring/2 +  tang_bins/2;
        else u = -abs(crystal1 - crystal2) +   detectors_per_ring/2 +  tang_bins/2;

        if ( u >=  tang_bins || u < 0 )
            continue;
        if (u%2 == 0)
        {
            zi = ( detectors_per_ring/2 - (crystal1 - crystal2) - 1)/2;
            if (zi >=   detectors_per_ring/4) zi = zi -  detectors_per_ring/2 + 1;
            if (zi <= - detectors_per_ring/4) zi = zi +  detectors_per_ring/2 - 1;
        }
        else
        {
            zi = ( detectors_per_ring/2 - (crystal1 - crystal2))/2;
            if (zi >=   detectors_per_ring/4) zi = zi -  detectors_per_ring/2;
            if (zi <= - detectors_per_ring/4) zi = zi +  detectors_per_ring/2;
        }

        c1 = crystal1 + zi;
        c2 = crystal2 - zi;
        if (c1 >=  detectors_per_ring) c1 = c1 -  detectors_per_ring;
        if (c1 < 0)      c1 = c1 +  detectors_per_ring;
        if (c2 >=  detectors_per_ring) c2 = c2 -  detectors_per_ring;
        if (c2 < 0)      c2 = c2 +  detectors_per_ring;

        if (c1 < c2) flip = 0;
        else         flip = 1;

        if (flip)
        {
            swap  = ring1;
            ring1 = ring2;
            ring2 = swap;
        }

        // Update the different arrays...
        //-------------------------------
        //***ALL EVENTS
        if(!load_test)//no out put now
        {
        float_michelogram[ring2][ring1][phi][u] += static_cast<float>(1.);}

        recorded_counts++;

        average_depth_of_interaction += globalPosX1;
    }
        average_depth_of_interaction /= recorded_counts;
        average_depth_of_interaction -=  inner_diameter/static_cast<float>(2.0);

}
void rootGen::showFirstaa(){
    FILE *fpa=fopen("myOutFile.csv","w");
    FILE *fpc=fopen("myOutFile50.csv","w");
    FILE *fpb=fopen("myOutFile.dat","wb");
    fwrite(float_michelogram,sizeof (float),(480*480),fpb);
    for (int k=0;k<100;k++){
    for(int i=0;i<480;i++){
        for(int j=0;j<480;j++)
        {
            fprintf(fpa,"%d ",static_cast<int>(float_michelogram[k][k][i][j]));
            if(k==50) fprintf(fpc,"%d ",static_cast<int>(float_michelogram[k][k][i][j]));
            //std::cout<<float_michelogram[50][50][i][j]<<" ";
        }
        //std::cout<<std::endl;
    }
    }
    fclose(fpa);
    fclose(fpb);
}

void rootGen::saveMichelogram()
{
    long segment_number = 0 ;
    long ring1=number_of_rings, ring2=number_of_rings;
    std::string output_data_name = output_template_name+".s";
    std::remove(output_data_name.c_str());
    FILE *output_michelogram_file=std::fopen(output_data_name.c_str(),"wb");

    if(output_michelogram_file!=nullptr) std::cout<<"FileOpened"<<std::endl;

    //save_as_projections = ui->checkBox_4->isChecked();
    save_as_projections=true;

    if (!save_as_projections)
    {

        if (data_size=="float")
        {
            unsigned long size=maximum_ring_difference*maximum_ring_difference;
            size*=( detectors_per_ring/2);
            size*= tang_bins;
            size*=sizeof(float);
            fwrite(float_michelogram,sizeof(float),(size), output_michelogram_file);
            std::cout<<"fwriteOK,Length="<<size<<std::endl;
        }

        else if (data_size=="integer")
        {
            fwrite(int_michelogram,sizeof(int),(selected_ring_difference*selected_ring_difference*( detectors_per_ring/2)* tang_bins), output_michelogram_file);
        }

        else if (data_size=="short integer")
        {
            fwrite(short_michelogram,sizeof(short),(selected_ring_difference*selected_ring_difference*( detectors_per_ring/2)* tang_bins), output_michelogram_file);
        }
    }
    else
    {
        for (qint16 i = 0 ; i < 2*maximum_ring_difference + 1 ; i++)
        {

            if (i <= maximum_ring_difference)
                segment_number =  number_of_rings - maximum_ring_difference + i;
            else
                segment_number =  number_of_rings + maximum_ring_difference - i;
            std::cout<<"segment_num="<<segment_number<<"\t";
            float *Proj;
            if(segment_number==0) continue;
            for (qint16 j = 0 ; j <  detectors_per_ring/2 ; j++)
            {

                if (data_size=="float")
                {
                    Proj=static_cast<float*>(malloc((sizeof (float)*segment_number*tang_bins)));

                    for (int k = 0 ; k < segment_number ; k++)
                    {
                        if (i <= maximum_ring_difference)
                            ring1 = k;
                        else
                            ring2 = k;

                        if (i <= maximum_ring_difference)
                            ring2 = ring1 + maximum_ring_difference - i;
                        else
                            ring1 = ring2 - maximum_ring_difference + i;

                        for (int l = 0 ; l <  tang_bins ; l++)
                            Proj[k*segment_number+l] = float_michelogram[ring2][ring1][j][l];
                    }
                    fwrite(Proj,sizeof(float),(segment_number* tang_bins), output_michelogram_file);
                    free(Proj);
                }

                else if (data_size=="integer")
                {

                    qint32 *Proj=static_cast<qint32*>(malloc((sizeof (qint32)*abs(segment_number)*abs(tang_bins))));
                    for (int k = 0 ; k < segment_number ; k++)
                    {
                        if (i <= selected_ring_difference)
                            ring1 = k;
                        else
                            ring2 = k;

                        if (i <= selected_ring_difference)
                            ring2 = ring1 + selected_ring_difference - i;
                        else
                            ring1 = ring2 - selected_ring_difference + i;

                        for (int l = 0 ; l <  tang_bins ; l++)
                            Proj[k*segment_number+l] = float_michelogram[ring2][ring1][j][l];
                    }
                    fwrite(Proj,sizeof(qint32),(segment_number* tang_bins), output_michelogram_file);
                }

                else if (data_size=="short integer")
                {
                    qint16 **Proj=static_cast<qint16**>(malloc((sizeof (qint16)*abs(segment_number)*abs(tang_bins))));

                    for (int k = 0 ; k < segment_number ; k++)
                    {
                        if (i <= selected_ring_difference)
                            ring1 = k;
                        else
                            ring2 = k;

                        if (i <= selected_ring_difference)
                            ring2 = ring1 + selected_ring_difference - i;
                        else
                            ring1 = ring2 - selected_ring_difference + i;

                        for (int l = 0 ; l <  tang_bins ; l++)
                            Proj[k][l] = float_michelogram[ring2][ring1][j][l];
                    }
                    fwrite(Proj,sizeof(qint16),(segment_number* tang_bins), output_michelogram_file);
                }
            }
        }
    }

    fclose(output_michelogram_file);

}
void rootGen::processMacData(){

    float vis_viewer_set_viewpoint_Theta=60,vis_viewer_set_viewpoint_Phi=60;
    float gate_world_geo_xLength=400;//cm
    float gate_world_geo_yLength=400;//cm
    float gate_world_geo_zLength=400;//cm
    float gate_PET_geo_Rmax=52;//cm
    float gate_PET_geo_Rmin=39.9f;//cm
    float gate_PET_geo_Height=40.2f;//cm
    float gate_head_geo_xLength=8;//cm
    float gate_head_geo_yLength=32;//cm
    float gate_head_geo_zLength=40;//cm
    float gate_block_geo_xLength=30;//mm
    float gate_block_geo_yLength=16;//mm
    float gate_block_geo_zLength=20;//mm
    float gate_crystal_geo_xLength=30;//mm
    float gate_crystal_geo_yLength=3;//mm
    float gate_crystal_geo_zLength=3.8f;//mm
    float gate_LSO_geo_xLength=15;//mm
    float gate_LSO_geo_yLength=3;//mm
    float gate_LSO_geo_zLength=3.8f;//mm
    float gate_BGO_geo_xLength=15;//mm
    float gate_BGO_geo_yLength=3;//mm
    float gate_BGO_geo_zLength=3.8f;//mm
    int gate_crystal_array_x=1;
    int gate_crystal_array_y=5;
    int gate_crystal_array_z=5;
    int gate_block_array_x=1;
    int gate_block_array_y=20;
    int gate_block_array_z=20;
    int gate_ring_count=4;
    //user defined params

    has_block=false;
    has_module=false;
    has_rsector=true;
    inner_ring_param1=2;//crystal_array.at(6)
    inner_ring_param2=1;inner_ring_param3=1;inner_ring_param4=1;inner_ring_param5=31.55;
    //crystal_array[3],block_array[3],module_array.at(3),rsector_array.at(1)
    number_of_detector_per_ring=960;
    //number_of_detector_per_ring =crystal_array.at(1) * block_array.at(1) * module_array.at(1) * rsector_array.at(0);
    maximum_ring_difference=number_of_rings-1;
    minimum_ring_difference=0;
    blocks_xy=1;
    blocks_z=1;
    crystals_z=20;
    crystals_xy=20;
    modules_z=5;
    modules_xy=1;
    number_of_rings=100;
    low_energy_u=450;
    hi_energy_u=650;
    view_offset_degrees=0;
    image_scaling_factor=0;
    time_frames=0;
    //number_of_detector_per_ring=gate_crystal_array_x*gate_block_array_x;

    //number_of_detector_per_ring =setRepeatNumberY *1*1*setRepeatNumber
    //param process
    detectors_per_ring=number_of_detector_per_ring;
    number_of_bins=number_of_detector_per_ring/2;
    tang_bins=number_of_bins;

    selected_ring_difference=maximum_ring_difference - minimum_ring_difference +1 ;
    max_ring_diff =number_of_rings;

    suggested_offset=number_of_detector_per_ring;
    gate_offset=static_cast<qint16> (number_of_detector_per_ring*0.75);


    //modification
    //int density_factor=2;
    //selected_ring_difference=1;
    //number_of_bins=number_of_detector_per_ring/2/density_factor;
    //tang_bins=number_of_bins/density_factor;
    //max_ring_diff=1;
}
int rootGen::loadROOTfiles()
{
    //bool ok;

    //qint32 max_num_of_events = 0;



    if(Coincidences!=nullptr) std::cout<<"exists";
    Coincidences->Add("file.root");

    Coincidences->SetBranchStatus("*",0);
    Coincidences->SetBranchAddress("comptonPhantom1",&comptonPhantom1);
    Coincidences->SetBranchAddress("comptonPhantom2",&comptonPhantom2);
    Coincidences->SetBranchAddress("sourceID1", &sourceID1);
    Coincidences->SetBranchAddress("sourceID2", &sourceID2);
    Coincidences->SetBranchAddress("eventID1",&eventID1);
    Coincidences->SetBranchAddress("eventID2",&eventID2);
    Coincidences->SetBranchAddress("energy1",&energy_1);
    Coincidences->SetBranchAddress("energy2",&energy_2);
    Coincidences->SetBranchAddress("submoduleID1",&submoduleID1);
    Coincidences->SetBranchAddress("submoduleID2",&submoduleID2);
    Coincidences->SetBranchAddress("crystalID1",&crystalID1);
    Coincidences->SetBranchAddress("crystalID2",&crystalID2);
    Coincidences->SetBranchAddress("moduleID1",&moduleID1);
    Coincidences->SetBranchAddress("moduleID2",&moduleID2);
    Coincidences->SetBranchAddress("rsectorID1",&rsectorID1);
    Coincidences->SetBranchAddress("rsectorID2",&rsectorID2);
    Coincidences->SetBranchAddress("globalPosX1",&globalPosX1);
    Coincidences->SetBranchAddress("globalPosX2",&globalPosX2);
    Coincidences->SetBranchAddress("sourcePosX1",&sourcePosX1);
    Coincidences->SetBranchAddress("sourcePosY1",&sourcePosY1);
    Coincidences->SetBranchAddress("sourcePosZ1",&sourcePosZ1);
    Coincidences->SetBranchAddress("sourcePosX2",&sourcePosX2);
    Coincidences->SetBranchAddress("sourcePosY2",&sourcePosY2);
    Coincidences->SetBranchAddress("sourcePosZ2",&sourcePosZ2);
    //combineROOTFiles();
    nentries = (double)(Coincidences->GetEntries());
    std::cout<<"nentries="<<nentries<<std::endl;
    //createROOTMichelogram();

    //find_maximums();
    return 1;


}

#ifdef debug
    void rootGen::debugWork(){
    /*Coincidences->compton1=&compton1;
    Coincidences->compton2=&compton2;
    Coincidences->runID=&runID;
    Coincidences->eventID1=&eventID1;
    Coincidences->eventID2=&eventID2;
    Coincidences->crystalID1=&crystalID1;
    Coincidences->crystalID2=&crystalID2;
    Coincidences.submoduleID1=&submoduleID1;
    Coincidences.submoduleID2=&submoduleID2;
    Coincidences.moduleID1=&moduleID1;
    Coincidences.moduleID2=&moduleID2;
    Coincidences.rsectorID1=&rsectorID1;
    Coincidences.rsectorID2=&rsectorID2;
    Coincidences.rotAngleID1=&rotAngleID1;
    Coincidences.rotAngleID2=&rotAngleID2;
    Coincidences.energy_1=&energy_1;
    Coincidences.energy_2=&energy_2;
    Coincidences.globalPosX1=&globalPosX1;
    Coincidences.globalPosX2=&globalPosX2;
    Coincidences.comptonPhantom1=&comptonPhantom1;
    Coincidences.comptonPhantom2=&comptonPhantom2;
    Coincidences.sourceID1=&sourceID1;
    Coincidences.sourceID2=&sourceID2;
    Coincidences.sourcePosX1=&sourcePosX1;
    Coincidences.sourcePosY1=&sourcePosY1;
    Coincidences.sourcePosZ1=&sourcePosZ1;
    Coincidences.sourcePosX2=&sourcePosX2;
    Coincidences.sourcePosY2=&sourcePosY2;
    Coincidences.sourcePosZ2=&sourcePosZ2;
    printf("linkOK\n");*/
    }
#endif



