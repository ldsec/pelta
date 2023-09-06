mkdir threads_experiment

go run main.go --n 10 --threads 2 --rel keygen > threads_experiment/t2.txt
go run main.go --n 10 --threads 4 --rel keygen > threads_experiment/t4.txt
go run main.go --n 10 --threads 8 --rel keygen > threads_experiment/t8.txt
go run main.go --n 10 --threads 16 --rel keygen > threads_experiment/t16.txt
go run main.go --n 10 --threads 32 --rel keygen > threads_experiment/t32.txt
go run main.go --n 10 --threads 48 --rel keygen > threads_experiment/t48.txt
go run main.go --n 10 --threads 64 --rel keygen > threads_experiment/t64.txt
go run main.go --n 10 --threads 96 --rel keygen > threads_experiment/t96.txt
