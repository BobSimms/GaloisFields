import sys
import random

def gcdex(a, b):
    x, y, r, xx, yy, rr = 1, 0, a, 0, 1, b
    while rr !=0:
        q = r //rr
        x, y, r, xx, yy, rr = xx, yy, rr, x -q *xx, y -q *yy, r -q *rr
    return x, y, r


class IntegersModP(object):
    def __init__(self, p):
        self.p = p
        self.order = self.p
        self.zero = 0
        self.one = 1
        self.elementDict = dict([(e, True) for e in range(self.p)])
        self.elements = tuple(range(self.p))
    def mul(self, a, b):
        return a *b %self.p
    def add(self, a, b):
        return (a +b) %self.p
    def sub(self, a, b):
        return (a -b) %self.p
    def div(self, a, b):
        s, t, gcd = gcdex(b, self.p)
        return (s *a) %self.p
    def __unicode__(self):
        return 'F' +str(self.p)
    def __str__(self):
        return unicode(self)
    def __contains__(self, x):
        return x in self.elementDict
    def __getitem__(self, i):
        return self.elements[i]
    def __len__(self):
        return self.order


class GF(object):
    """Provides arithmetic for a field where the elements are polynomial remainders for polynomials having coefficients from any specified field.  The number of elements is a prime-power"""

    def __init__(self, cfield, dimension, symbol='x', modPoly='random'):
        self.cf = cfield
        self.dimension = dimension # number of terms or coordinates
        self.order = self.cf.order **self.dimension
        self.zero = (self.cf.zero,) *self.dimension
        self.one = (self.cf.one,) +(self.cf.zero, ) *(self.dimension -1)
        self.modPoly = tuple(modPoly)
        self.symbol = symbol

        while True:
            if modPoly =='random':
                self.modPoly = tuple(random.choice(self.cf.elements) for i in range(self.dimension)) +(self.cf.one,)
            assert(len(self.modPoly) == self.dimension +1)
            try:
                 self.build()
            except RuntimeError:
                if modPoly =='random':
                    continue
                else:
                    raise
            break

    def build(self):
        self.log = {}
        self.alog = [self.cf.one] *(self.order -1)
        self.elementDict = {self.zero:True}
        val = self.one
        for i in range(0, self.order -1):
            if val in self:
                raise RuntimeError(("ERROR in GF log/alog tables, 2^%d = %s is a repeat value.\n"
                +"   The characteristic polynomial %s must be invalid") % (i +1, val, str(self.modPoly)))
            self.elementDict[val] = True
            self.log[val] = i
            self.alog[i] = val
            val = self.normalize((self.cf.zero, ) +val) # mult by (0, 1) = x^1
        self.elements = tuple([self.zero] +self.alog)

    def __unicode__(self):
        if type(self.cf).__name__ =='IntegersModP':
            fmt = '%s[%s]/%s'
        else:
            fmt = '(%s)[%s]/%s'
        if type(self.cf).__name__ =='IntegersModP':
            return fmt % (unicode(self.cf), self.symbol, unicode(self.modPoly))
        else:
            return fmt % (unicode(self.cf), self.symbol, '(' +', '.join([(self.cf.symbol +'^' +str(self.cf.log[c]) if c!=self.cf.zero else '0') for c in self.modPoly]) +')')
    def __str__(self):
        return unicode(self)
    def __contains__(self, x):
        return x in self.elementDict
    def __getitem__(self, i):
        return self.elements[i]
    def __len__(self):
        return self.order

    def normalize(self, x):
        if len(x) >self.dimension:
            m = x[self.dimension]
            if m !=self.cf.zero:
                temp = tuple(self.cf.sub(x[i], self.cf.mul(self.modPoly[i], m)) for i in range(self.dimension))
            else:
                temp = x[:self.dimension]
        else:
            temp = x
        return temp

    def mul(self, x, y):
        if x==self.zero or y ==self.zero:
            return self.zero
        else:
            return self.alog[(self.log[x] +self.log[y]) %(self.order -1)]

    def pow(self, x, n):
        if x ==self.zero:
            return self.zero
        else:
            return self.alog[(self.log[x] *n) %(self.order -1)]

    def div(self, x, y):
        if y ==self.zero:
            raise ZeroDivisionError
        elif x ==self.zero:
            return self.zero
        else:
            return self.alog[(self.log[x] -self.log[y]) %(self.order -1)]

    def add(self, x, y):
        return tuple(self.cf.add(a, b) for a, b in zip(x, y))

    def sub(self, x, y):
        return tuple(self.cf.sub(a, b) for a, b in zip(x, y))

    def info(self, style='quotient'):
        # style = quotient poly coefficients or extension powers, thereof
        print(unicode(self))
        cols = 2
        rows = ((self.order -1) +1) //2
        for row in range(rows):
            for col in range(cols):
                index = col *rows +row
                if index <= self.order -2:
                    if style =='extension' and type(self.cf).__name__ =='GF':
                        sys.stdout.write('    %s^%d = (%s)' % (self.symbol, index, ', '.join([(self.cf.symbol +'^' +str(self.cf.log[c]) if c!=self.cf.zero else '0') for c in self.alog[index]])))
                    elif style =='quotient':
                        sys.stdout.write('    %s^%d = %s' % (self.symbol, index, self.alog[index]))
            sys.stdout.write('\n')
        print('')
