def factorial(n):
    if n==1:
        return 1
    else:
        answer=n*factorial(n-1)
        return answer


number=int(input("Enter a number: "))
print(f"the factorial of {number} is {factorial(number)}")


