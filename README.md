# LattiCom
An experimental library for lattice-based commitments

## Installation
The only requirement is Golang. The dependencies are pulled automatically. The following shows how to run the code on an Ubuntu based distribution.
```shell
$ cd LattiCom
$ sudo apt install -y golang
$ go run main.go [PARAMETERS]
```

## Parameters
- `--help`: Displays information about the parameters.
- `--n`: Number of executions. Defaults to `1`.
- `--rlwe`: Log RLWE ring degree. This can be `10`, `11`, `12`, `13`, `14`, or `15`. Defaults to `13`.
- `--commt`: Log commitment ring degree. Defaults to `7`.
- `--levels`: Number of levels for the ring. This can be `1`, `2`, or `3`. Defaults to `1`.
- `--delta1`: Log2 delta1. Defaults to `25`.
- `--lambda`: Lambda security parameter. Defaults to `10`.
- `--kappa`: Kappa security parameter. Defaults to `9`.
- `--rel`: The relation to run. Can be `keygen`, `aggr_keygen`, `collective_dec`, `keyswitch`, `bootstrapping`, or `relinkeygen`. Defaults to `keygen`.

## Interpreting the output
The program runs the protocol `n` many times where `n` is the number of executions supplied through the program arguments.
The final output contains both the average and individual running times for `Setup`, `Prover`, and `Verifier` across the runs.
For every data point, we also report the timings of its subprocedures. All the presented numbers are in milliseconds.

## Example configurations
Here we give the run configuration to be able to run the experiments we have performed for our paper.

**Important:** Note that the files need to be regenerated for every ring. Please remove the cached files in the root
folder by invoking `rm *.test` after switching the ring parameters.
### Table 1: Key generation
```shell
$ go run main.go --n 10 --rel keygen
```
### Table 2: Other relations
```shell
$ go run main.go --n 10 --rel collective_dec
```
```shell
$ go run main.go --n 10 --rel keyswitch
```
```shell
$ go run main.go --n 10 --rel bootstrapping
```
### Table 3: Collective key generation
```shell
$ go run main.go --n 10 --rel aggr_keygen
```
### Table 4: Varying RLWE ring degree
```shell
$ go run main.go --n 10 --rlwe 11
```
```shell
$ go run main.go --n 10 --rlwe 12
```
```shell
$ go run main.go --n 10 --rlwe 13
```
```shell
$ go run main.go --n 10 --rlwe 14
```
```shell
$ go run main.go --n 10 --rlwe 15
```
### Table 5: Varying levels
```shell
$ go run main.go --n 10 --levels 1
```
```shell
$ go run main.go --n 10 --levels 2
```
```shell
$ go run main.go --n 10 --levels 3
```
### Table 6: Varying commitment ring degree
```shell
$ go run main.go --n 10 --commt 7 --lambda 10 --kappa 9 --delta1 25
```
```shell
$ go run main.go --n 10 --commt 10 --lambda 1 --kappa 2 --delta1 29
```
```shell
$ go run main.go --n 10 --commt 13 --lambda 1 --kappa 1 --delta1 32
```
