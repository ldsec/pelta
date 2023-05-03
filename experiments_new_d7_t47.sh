
mkdir tbl1 tbl2 tbl3 tbl4 tbl5 tbl6

go run main.go --n 10 --threads 47 --commt 7 --lambda 17 --kappa 8 --rel keygen > tbl1/keygen.txt
go run main.go --n 10 --threads 47 --commt 7 --lambda 17 --kappa 8 --rel collective_dec > tbl2/collective_dec.txt
go run main.go --n 10 --threads 47 --commt 7 --lambda 17 --kappa 8 --rel keyswitch > tbl2/keyswitch.txt
go run main.go --n 10 --threads 47 --commt 7 --lambda 17 --kappa 8 --rel bootstrapping > tbl2/bootstrapping.txt
go run main.go --n 10 --threads 47 --commt 7 --lambda 17 --kappa 8 --rel relinkeygen > tbl2/relinkeygen.txt
go run main.go --n 10 --threads 47 --commt 7 --lambda 17 --kappa 8 --rel aggr_keygen > tbl3/aggr_keygen.txt
rm *.test
go run main.go --n 10 --threads 47 --commt 7 --lambda 17 --kappa 8 --rlwe 11 > tbl4/keygen11.txt
rm *.test
go run main.go --n 10 --threads 47 --commt 7 --lambda 17 --kappa 8 --rlwe 12 > tbl4/keygen12.txt
rm *.test
go run main.go --n 10 --threads 47 --commt 7 --lambda 17 --kappa 8 --rlwe 13 > tbl4/keygen13.txt
rm *.test
go run main.go --n 10 --threads 47 --commt 7 --lambda 17 --kappa 8 --rlwe 14 > tbl4/keygen14.txt
rm *.test
go run main.go --n 10 --threads 47 --commt 7 --lambda 17 --kappa 8 --rlwe 15 > tbl4/keygen15.txt
rm *.test
go run main.go --n 10 --threads 47 --commt 7 --lambda 17 --kappa 8 --levels 1 > tbl5/keygen1.txt
rm *.test
go run main.go --n 10 --threads 47 --commt 7 --lambda 17 --kappa 8 --levels 2 > tbl5/keygen2.txt
rm *.test
go run main.go --n 10 --threads 47 --commt 7 --lambda 17 --kappa 8 --levels 3 > tbl5/keygen3.txt
rm *.test
go run main.go --n 10 --threads 47 --commt 7 --lambda 17 --kappa 8 --delta1 25 > tbl6/keygen7.txt
rm *.test
go run main.go --n 10 --threads 47 --commt 10 --lambda 3 --kappa 2 --delta1 29 > tbl6/keygen10.txt
rm *.test
go run main.go --n 10 --threads 47 --commt 13 --lambda 1 --kappa 1 --delta1 32 > tbl6/keygen13.txt
rm *.test

