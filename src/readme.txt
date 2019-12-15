to run sequential:

$ mpicc -o sequential.exe knnring_sequential.c tester.c -lm -lopenblas && ./sequential.exe

to run mpi synchronous:

$ mpicc -o sync.exe knnring_synchronous.c tester_mpi.c -lm -lopenblas && mpirun -np 4 ./sync.exe

to run mpi asynchronous:

$ mpicc -o async.exe knnring_asynchronous.c tester_mpi.c -lm -lopenblas && mpirun -np 4 ./async.exe
