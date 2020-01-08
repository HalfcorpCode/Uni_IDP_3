# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 01:42:18 2020

@author: The Halflife390
"""

for i in range(101):
    if i%3 == 0:
        if i%5 == 0:
            print("FizzBuss")
        else:
            print("Fizz")
    elif i%5 == 0:
        print("Buzz")
    else:
        print(str(i))
        