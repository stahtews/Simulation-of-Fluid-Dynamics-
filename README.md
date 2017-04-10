# Simulation-of-Fluid-Dynamics-
Parallelizing Wave Diffusion with MPI and OpenMP

The focus is to parallelize the serial 2D Wave program which diffuses along x and y axis
using Hybrid of MPI and OpenMP strategies . The height of Wave at any point in time t is calculated
using Schroedinger’s wave dissemination algorithm by using Z[t-1], Z[t-2] and neighboring point
values.

Serial Code Strategy : As the current value of the matrix depends on only t-1,t-2 values according to
the formula , we use just 3 matrices to represent data and rotated in next iteration. The value is
calculated for each time slot and intermediate Wave height is printed at intervals as specified by the
user.

The Parallelization strategy mainly depends on the partitioning of the entire Wave representation
into N stripes either along x axis or y axis .where N is the Value of number of nodes or processors
across the network which compute the values of respective stripe.

In MPI parallelization strategy, a ring of nodes neatly runs the local copies of the program. Special
instructions for all the nodes, including Master node (rank =0) are specified by defining group of “if
blocks” based on the ranks of nodes. Each node identified by my_rank will execute those
instructions to send, receive or process the data. In our current program, a separate strip of data
along x axis or y axis is assigned to each node. As mentioned earlier the value at Z[t] depends on Z[t-
1], Z[t-2] and neighboring points values. As each of the processor is already exposed to the values of
Z[t-1], z[t-2] from the beginning and hence requires only the start, end indexes of each stripe and
border values. The Indexes are dynamically calculated using my_rank, stripe, size etc. The
neighboring stripes are sent over the MPI_COMM_WORLD communication channel to next and
previous processors. These values are used for Schrodinger’s calculation. The intermediate results
are printed on demand by master node by gathering data using MPI_Send AND MPI_Recv, based on
the interval values specifies by the user. The Z[t] value at current time t, becomes Z[t-1] in the next
iteration and this process continues till the end of Max_time . This parallelization technique involves
communication over MPI_COMM_WORLD and hence the communicated data size and
synchronization of data is kept to minimum (only indexes and border stipes are sent and very less
synchronization) as it causes overhead.

OpenMP Parallelization strategy depends on partitioning of for loops based on the directives by
operating system. Threads are created for each processor, and the load is distributed equally among all
threads as the code is embarrassingly parallel because of lack of need for synchronization between 
iterations. These threads follow shared memory strategy and hence the private data or stack memory
usage is kept to minimum for each instance.

Combining these two techniques results in N * n (N – distributed processors, n – shared memory
processors) instances, parallel working on the local copy of the program. For example 4 processors and 4
threads in each processor creates 16 instances of individually working processors, which greatly impacts
the performance of the program.
