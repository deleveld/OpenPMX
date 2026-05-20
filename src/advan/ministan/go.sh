#!/bin/sh

rm ./main
gcc main_de.c -I. -lm -o main && ./main
