#!/usr/bin/env python

import sys

print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)

def sum_digits(input): #input is an integer
    result = 0

    input_list = []
    input_string = str(input)
    for letter in input_string:
        input_list.append(int(letter))

    for element in input_list:
        result += element

    return result




input1 = int(sys.argv[1])
if(len(sys.argv) > 2):
    input2 = int(sys.argv[2])
    printval = input1
    while(printval <= input2):
        print sum_digits(printval)
        printval += input1
else:
    printval = input1
    
    while(printval <= input1**2):
        print sum_digits(printval)
        printval += input1 #calculate next integer multiple
