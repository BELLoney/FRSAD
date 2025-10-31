clc;
clear all;
format short

load Example.mat

Dataori=trandata;


delta=1;
sigma=0.1;
out_factors=FRSAD(trandata,delta,sigma)