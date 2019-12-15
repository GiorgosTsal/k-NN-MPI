# k-NN ring MPI validators

> This is an automated validation case based on testers our teachers gave us in order to validate our results. To check the validity of our code you
can download the validation folder with all the essentials and run the following commands.

## Run the testers

To run the automated test-suite of some kind, do this for each version you want to evaluate.

* To run sequential version
```sh
make test_sequential
```
* To run distributed MPI ring synchronous version
```sh
make test_mpi_sync
```
* To run distributed MPI ring asynchronous version
```sh
make test_mpi_async
```



## Meta

Giorgos Tsalidis – gtsalidis@ece.auth.gr
Giannis Loias – loiasioann@ece.auth.gr

Distributed under the GNU General Public License v3.0. See ``LICENSE`` for more information.
