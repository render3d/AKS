/*
    A C++ source code by Yves Gallot that uses the NTL library for polynomial multiplication (FFT, combined with the Chinese Remainder Theorem):
        http://yves.gallot.pagesperso-orange.fr/src/aks_ntl.html
*/

#include "NTL/ZZ_pX.h"
#include <math.h>

NTL_CLIENT

#ifdef _M_IX86

#define umulrem(z, x, y, m)	\
__asm	mov		eax, x	\
__asm	mul		y	\
__asm	div		m	\
__asm	mov		z, edx

#else

#ifdef _MSC_VER
typedef unsigned __int64	Tu64;
#else
typedef unsigned long long	Tu64;
#endif

#define umulrem(z, x, y, m)	\
	{	\
	z = (unsigned int)(x * (Tu64)y % m);	\
	}

#endif

static bool IsPrime(unsigned int n)
{
	if (n < 2) return false;
	if (n < 4) return true;
	if (n % 2 == 0) return false;

	const unsigned int iMax = (int)sqrt(n) + 1;
	unsigned int i;
	for (i = 3; i <= iMax; i += 2)
		if (n % i == 0)
			return false;

	return true;
}

static unsigned int LargestPrimeFactor(unsigned int n)
{
	if (n < 2) return 1;

	unsigned int r = n, p;
	if (r % 2 == 0)
	{
		p = 2;
		do { r /= 2; } while (r % 2 == 0);
	}
	unsigned int i;
	for (i = 3; i <= r; i += 2)
	{
		if (r % i == 0)
		{
			p = i;
			do { r /= i; } while (r % i == 0);
		}
	}
	return p;
}

static unsigned int Powm(unsigned int n, unsigned int e, unsigned int m)
{
	unsigned int r = 1;
	unsigned int t = n % m;
	unsigned int i;
	for (i = e; i != 0; i /= 2)
	{
		if (i % 2 != 0)
		{
			umulrem(r, r, t, m);
		}
		umulrem(t, t, t, m);
	}
	return r;
}

int main(int argc, char * argv[])
{
// AKS
// n=11701, r=11699, q=5849, s=2923
// n=1000000007, r=57287, q=28643, s=14311
//
// Bernstein
// n=349, r=347, q=173, s=140
// n=1000000007, r=3623, q=1811, s=1785

	unsigned int n0;
	cout << "n ? ";
	cin >> n0;

	unsigned int n;
	for (n = n0; true; n += 2)
	{
		bool b = false;
		unsigned int q, r, s;
		for (r = 3; r < n; r += 2)
		{
			if (IsPrime(r))
			{
				const unsigned int m = n % r;
				if (m == 0) break;
				q = LargestPrimeFactor(r - 1);
				if (Powm(m, (r - 1) / q, r) <= 1) continue;
/*
				// AKS
				if (q >= 4 * sqrt(r) * n.Log()/log(2))
				{
					s = (unsigned int)(2 * sqrt(r) * n.Log()/log(2));
					b = true;
					break;
				}
*/
				// Bernstein
				const double cMin = 2 * floor(sqrt(r)) * log(n);
				double c = 0;
				for (s = 1; s <= q; s++)
				{
					c += log(q + s - 1) - log(s);
					if (c >= cMin)
					{
						b = true;
						break;
					}
				}
				if (b) break;

			}
		}
		if (b)
		{
			cout << "n=" << n << ", r=" << r << ", q=" << q << ", s=" << s << "\r" << flush;

			bool b = true;

			ZZ_p::init(to_ZZ(n));
			ZZ_pX f(r, 1); f -= 1;
			const ZZ_pXModulus pf(f);

			ZZ_pX rPoly(1, 1);				// x
			PowerMod(rPoly, rPoly, n, pf);	// x^n

			unsigned int a;
			for (a = 1; a <= s; ++a)
			{
				cout << "n=" << n << ", r=" << r << ", q=" << q << ", s=" << s << " " << a << "\r" << flush;;

				ZZ_pX lPoly(1, 1);	// x
				lPoly -= a;			// x - a
				PowerMod(lPoly, lPoly, n, pf);	// (x - a)^n
				lPoly += a;			// (x - a)^n + a
				if (lPoly != rPoly)
				{
					b = false;
					break;
				}
			}

			if (b)
				cout << "n=" << n << ", r=" << r << ", q=" << q << ", s=" << s << ", power of a prime." << endl;
			else
				cout << "n=" << n << ", r=" << r << ", q=" << q << ", s=" << s << ", composite." << endl;
		}
	}

	return 0;
}