---
title: "AOC: Day 1"
author: "FA"
date: 
tags:
categories:
---

## Problem:
Find the pair of numbers from a input that sum upto a certain number


```python
# the test data is sorted in ascending order
test = [1721,979,366,299,675,1456]
test = sorted(test)
test
```




    [299, 366, 675, 979, 1456, 1721]



## Solution Approach:
Algorithm:
1. Sort the provided data.  
2. two indices pointing to the low and high end of the sorted input
3. sum the values indicating the low and high indices
4. if the sum is less than the total, then increase the low indices
5. else if the sum is higher that the total, then decreaset the high indices
6. repeat step 3 to 5 until the solution is reached.


```python
def find_number(A, sum):
    """this function finds a number (sum) that is a summation of two numbers in the list (A)"""
    # sort the list in ascending order
    A.sort()
 
    # maintain two indices pointing to end-points of the list
    (low, high) = (0, len(A) - 1)
 
    # reduce search space A[low..high] at each iteration of the loop
 
    # loop till low is less than high
    while low < high:
 
        if A[low] + A[high] == sum:        # sum found
            print("Pair found")
            print(A[low], A[high])
            return
 
        # increment low index if total is less than the desired sum
        # decrement high index is total is more than the sum
        if A[low] + A[high] < sum:
            low = low + 1
        else:
            high = high - 1
  
    
```


```python
# test if the function works or not
find_number([1721,979,366,299,675,1456], 2020)
```

    Pair found
    299 1721
    


```python
# read the data

my_file = open("input.txt", "r")
content = my_file.read()
content_list = content.split("\n")[:-1]
my_file.close()

data = [int(i) for i in content_list]
print(data)
```

    [1713, 1281, 1185, 1501, 1462, 1752, 1363, 1799, 1071, 1446, 1685, 1706, 1726, 1567, 1867, 1376, 1445, 1971, 1429, 1749, 438, 1291, 1261, 1585, 1859, 1835, 1630, 1975, 1467, 1829, 1669, 1638, 1961, 1719, 1238, 1751, 1514, 1744, 1547, 1677, 1811, 1820, 1371, 740, 1925, 1803, 1753, 1208, 1772, 1642, 1140, 1838, 1444, 1321, 1556, 1635, 1687, 688, 1650, 1580, 1290, 1812, 1814, 1384, 1426, 1374, 1973, 1791, 1643, 1846, 1676, 1724, 1810, 1911, 1765, 945, 1357, 1919, 1994, 1697, 1632, 1449, 1539, 1725, 1963, 1879, 1731, 1904, 1392, 1823, 1420, 1504, 204, 1661, 1575, 1401, 1806, 1417, 1965, 1960, 1990, 1409, 1649, 1566, 1957, 514, 1464, 1352, 1841, 1601, 1473, 1309, 1421, 1190, 1582, 1825, 655, 1666, 1878, 1891, 1579, 1176, 1557, 1910, 1747, 1388, 1493, 1372, 1522, 1515, 1745, 1494, 1763, 1147, 1364, 1469, 1165, 1901, 1368, 1234, 1308, 1416, 1678, 1541, 1509, 1427, 1223, 1496, 1600, 1383, 1295, 1415, 1890, 1694, 1793, 1529, 1984, 1576, 1244, 1348, 1085, 1770, 1358, 1611, 1159, 1964, 1647, 818, 1246, 1458, 1936, 1370, 1659, 1923, 1619, 1604, 1354, 1118, 1657, 1945, 1898, 1948, 798, 769, 1689, 1821, 1979, 1460, 1832, 1596, 1679, 1818, 1815, 1977, 1634, 1828, 1386, 1284, 1569, 1970]
    


```python
find_number(data, 2020)
```

    Pair found
    438 1582
    


```python
for i in data:
    print(f"i:{i}")
    find_number(data, 2020-i)
```

    i:204
    i:438
    i:514
    Pair found
    688 818
    i:655
    i:688
    Pair found
    514 818
    i:740
    i:769
    i:798
    i:818
    Pair found
    514 688
    i:945
    i:1071
    i:1085
    i:1118
    i:1140
    i:1147
    i:1159
    i:1165
    i:1176
    i:1185
    i:1190
    i:1208
    i:1223
    i:1234
    i:1238
    i:1244
    i:1246
    i:1261
    i:1281
    i:1284
    i:1290
    i:1291
    i:1295
    i:1308
    i:1309
    i:1321
    i:1348
    i:1352
    i:1354
    i:1357
    i:1358
    i:1363
    i:1364
    i:1368
    i:1370
    i:1371
    i:1372
    i:1374
    i:1376
    i:1383
    i:1384
    i:1386
    i:1388
    i:1392
    i:1401
    i:1409
    i:1415
    i:1416
    i:1417
    i:1420
    i:1421
    i:1426
    i:1427
    i:1429
    i:1444
    i:1445
    i:1446
    i:1449
    i:1458
    i:1460
    i:1462
    i:1464
    i:1467
    i:1469
    i:1473
    i:1493
    i:1494
    i:1496
    i:1501
    i:1504
    i:1509
    i:1514
    i:1515
    i:1522
    i:1529
    i:1539
    i:1541
    i:1547
    i:1556
    i:1557
    i:1566
    i:1567
    i:1569
    i:1575
    i:1576
    i:1579
    i:1580
    i:1582
    i:1585
    i:1596
    i:1600
    i:1601
    i:1604
    i:1611
    i:1619
    i:1630
    i:1632
    i:1634
    i:1635
    i:1638
    i:1642
    i:1643
    i:1647
    i:1649
    i:1650
    i:1657
    i:1659
    i:1661
    i:1666
    i:1669
    i:1676
    i:1677
    i:1678
    i:1679
    i:1685
    i:1687
    i:1689
    i:1694
    i:1697
    i:1706
    i:1713
    i:1719
    i:1724
    i:1725
    i:1726
    i:1731
    i:1744
    i:1745
    i:1747
    i:1749
    i:1751
    i:1752
    i:1753
    i:1763
    i:1765
    i:1770
    i:1772
    i:1791
    i:1793
    i:1799
    i:1803
    i:1806
    i:1810
    i:1811
    i:1812
    i:1814
    i:1815
    i:1818
    i:1820
    i:1821
    i:1823
    i:1825
    i:1828
    i:1829
    i:1832
    i:1835
    i:1838
    i:1841
    i:1846
    i:1859
    i:1867
    i:1878
    i:1879
    i:1890
    i:1891
    i:1898
    i:1901
    i:1904
    i:1910
    i:1911
    i:1919
    i:1923
    i:1925
    i:1936
    i:1945
    i:1948
    i:1957
    i:1960
    i:1961
    i:1963
    i:1964
    i:1965
    i:1970
    i:1971
    i:1973
    i:1975
    i:1977
    i:1979
    i:1984
    i:1990
    i:1994
    


```python

```
