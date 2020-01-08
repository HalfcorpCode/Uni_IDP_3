"""
FizzBuzz Program
Gianluca Cantone 2019
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
        