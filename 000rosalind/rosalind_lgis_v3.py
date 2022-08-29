# Dynamic programming Python implementation of LIS problem

# lis returns length of the longest increasing subsequence
# in arr of size n
def lis(arr):
    n = len(arr)
    # Declare the list (array) for LIS and initialize LIS
    # values for all indexes
    lis = [1]*n

    prev = [0]*n
    for i in range(0, n):
        prev[i] = i

    # Compute optimized LIS values in bottom up manner
    for i in range (1 , n):
        for j in range(0 , i):
             #if int(arr[j]) < int(arr[i]) 
             if int(arr[i]) < int(arr[j]) and lis[i] < lis[j] + 1 :
             #if int(arr[i]) > int(arr[j]) and lis[i] < lis[j] + 1 :
                lis[i] = lis[j]+1
                #lis[i] = max(lis[i], lis[j]+1)
                prev[i] = j
    #print(prev), the prev of location i 
    #n log n

    # Initialize maximum to 0 to get the maximum of all
    # LIS
    maximum = 0
    idx = 0

    # Pick maximum of all LIS values
    for i in range(n):
        if maximum < lis[i]:
            maximum = lis[i]
            idx = i
        #print(lis, maximum, idx)

    seq = [arr[idx]]
    while idx != prev[idx]:
        #print(idx, prev[idx], seq)
        idx = prev[idx]
        seq.append(arr[idx])
        #print(idx, prev[idx], seq)
        # 5 4 2 not 5 4 3 

    return (maximum, reversed(seq))
# end of lis function

# Driver program to test above function
#arr = "5 1 4 2 3".split()
ans = lis(arr)
print("The longest increase sequence is:\n", " ".join(str(x) for x in ans[1]),"\n")
#rra = arr[::-1]
#ans = lis(rra)
#print("The longest decrease sequence is:\n", " ".join(str(x) for x in ans[1]),"\n")