//#####################################################################
// Copyright 2010, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __Kernel_Serial_Base_Helper__
#define __Kernel_Serial_Base_Helper__

#include <Thread_Queueing/PTHREAD_QUEUE.h>

using namespace PhysBAM;
extern PTHREAD_QUEUE* pthread_queue;

namespace MT_Streaming_Kernels{

    template<class OP> class Kernel_Serial_Base_Helper;

    namespace{
        template<class OP>
            struct Kernel_Serial_Base_Thread_Helper:public PTHREAD_QUEUE::TASK
            {
                Kernel_Serial_Base_Helper<OP>* const obj;
                const int entry_start,entry_end;
            Kernel_Serial_Base_Thread_Helper(Kernel_Serial_Base_Helper<OP>* const obj_input,const int entry_start_input,const int entry_end_input)
        :obj(obj_input),entry_start(entry_start_input),entry_end(entry_end_input) {}
                void Run(){obj->Run_Index_Range(entry_start,entry_end);}
            };
    }
    

    template<class OP>
        class Kernel_Serial_Base_Helper
        {
            OP& op;
            const int number_of_entries;
            const int number_of_partitions;

        public:
            explicit Kernel_Serial_Base_Helper(OP& op_input,
                                             const int number_of_entries_input,
                                             const int number_of_partitions_input
                                        )
                :op(op_input),
                number_of_entries(number_of_entries_input),
                number_of_partitions(number_of_partitions_input)
                {}
    
            void Run()
            {
                Run_Index_Range(0,number_of_entries);
            }
  
//#####################################################################

            void Run_Parallel()
            { 
                int partition_offsets[number_of_partitions];
                int partition;
                partition_offsets[0]=0;
                for( int i=0, partition=0; i < number_of_entries && partition < number_of_partitions;
                     i+=number_of_entries / number_of_partitions, partition++ )
                    partition_offsets[partition] = i;

                for(partition=0;partition<number_of_partitions;partition++){
                    int entry_begin=partition_offsets[partition];
                    int entry_end=((partition<number_of_partitions-1)?partition_offsets[partition+1]:number_of_entries);
            
                    Kernel_Serial_Base_Thread_Helper<OP>* task=
                        new Kernel_Serial_Base_Thread_Helper<OP>(this,entry_begin,entry_end);
                    pthread_queue->Queue(task);
                }
                pthread_queue->Wait();
            }

            void Run_Parallel_OMP()
            { 
                Run_Index_Range_OMP(0,number_of_entries,number_of_partitions);
            }

//#####################################################################


            void Run_Index_Range(const int entry_start,const int entry_end)
            {
                for(int entry=entry_start;entry<entry_end;entry++)
                    op.Execute( entry ); 
            }

            void Run_Index_Range_OMP(const int entry_start,const int entry_end, const int number_of_partitions)
            {
                #pragma omp parallel for num_threads(number_of_partitions)
                for(int entry=entry_start;entry<entry_end;entry++)
                    op.Execute( entry ); 
            }

//#####################################################################
        };
}
#endif
