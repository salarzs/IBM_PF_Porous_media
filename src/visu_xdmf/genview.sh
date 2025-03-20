#!/bin/bash

gfortran gen_xdmf.f90 && ./a.out && rm -rf a.out
