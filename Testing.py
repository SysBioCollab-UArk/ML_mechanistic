def average_of_x_numbers(lst, x):
    result = []
    for i in range(0, len(lst)-x+1):
        avg = sum(lst[i:i + x]) / x
        result.append(avg)
    return result

numbers = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
x = 2
print(average_of_x_numbers(numbers, x))