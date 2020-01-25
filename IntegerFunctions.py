def gcdex(a, b):
    x, y, r, xx, yy, rr = 1, 0, a, 0, 1, b
    while rr !=0:
        q = r //rr
        x, y, r, xx, yy, rr = xx, yy, rr, x -q *xx, y -q *yy, r -q *rr
    return x, y, r

def phi(n):
    # Initialize result as n
    result = n

    # Consider all prime factors 
    # of n and subtract their 
    # multiples from result 
    p = 2
    while(p * p <=n): 
        # Check if p is a
        # prime factor.
        if (n % p ==0):
            # If yes, then
            # update n and result
            while (n % p ==0):
                n = n //p
            result -= result //p
        p += 1

    # If n has a prime factor
    # greater than sqrt(n)
    # (There can be at-most
    # one such prime factor)
    # i.e. n is prime, so phi(n) = n -(n //n) = n -1
    if (n > 1):
        result -= result //n
    return result

def primeFactors(n):
    factors = []
    if n in [2, 3]:
        factors.append(n)
        return factors
    d=2
    while(d *d <=n):
        while(n >1):            
            while n %d ==0:
                factors.append(d)
                n = n /d
            d += 1
    return factors


if __name__ =='__main__':
    print('phi(100) = %d' % phi(100))
    print('gcdex(15, 11) = %s' % str(gcdex(15, 11)))
    print(primeFactors(5 *365))
