// Octave module for the low-level induced velocity computations that 
// consume the largest fraction of time in the nonlinear lifting-line 
// code. Using this module speeds up the computation at least 6x. 
// 
// Compile into a .oct file using:
// 
// macOS:
// mkoctfile "-O3 -march=core-avx2 -mavx2 -ffp-contract=fast -DNDEBUG" cbiot.cpp
//
// or, if your processor does not have support for AVX2 instructions:
// 
// mkoctfile "-O3 -march=core2 -DNDEBUG" cbiot.cpp
//
// Windows/Linux:
// mkoctfile "-O3 -march=native" cbiot.cpp
//  
// (c) 2017 David Eller david@larosterna.com 

#include <octave/oct.h>
#include <limits>

#pragma STDC FP_CONTRACT ON

class PListView {
public:
  PListView(double *p, size_t n) : m_data(p), m_npts(n) {}
  double &x(size_t j) { return m_data[j]; }
  const double &x(size_t j) const { return m_data[j]; }
  double &y(size_t j) { return m_data[m_npts+j]; }
  const double &y(size_t j) const { return m_data[m_npts+j]; }
  double &z(size_t j) { return m_data[2*m_npts+j]; }
  const double &z(size_t j) const { return m_data[2*m_npts+j]; }
  void fill(double v) {
    std::fill(m_data, m_data+3*m_npts, v);
  }
private:
  double *m_data;
  size_t m_npts;
};

class PRef 
{
public:
  PRef(PListView list, size_t idx) : m_list(list), m_idx(idx) {}
  PRef &operator= (double v) {x() = v; y() = v; z() = v; return *this;}
  double &x() { return m_list.x(m_idx);}
  const double &x() const { return m_list.x(m_idx);}
  double &y() { return m_list.y(m_idx);}
  const double &y() const { return m_list.y(m_idx);}
  double &z() { return m_list.z(m_idx);}
  const double &z() const { return m_list.z(m_idx);}
  void cross(const PRef &a, const PRef &b) {
    x() = a.y()*b.z() - a.z()*b.y();
    y() = a.z()*b.x() - a.x()*b.z();
    z() = a.x()*b.y() - a.y()*b.x();
  }
  void sub(const PRef &a, const PRef &b) {
    x() = a.x() - b.x();
    y() = a.y() - b.y();
    z() = a.z() - b.z();
  }
  void scale(double a) {
    x() *= a;
    y() *= a;
    z() *= a;
  }
  void print(const char *pfx) const {
    octave_stdout << pfx << " [" << x() << ", " << y() << ", " << z() << "]\n";
  }
private:
  PListView m_list;
  size_t m_idx;
};

static inline double dot(const PRef &a, const PRef &b)
{
  return a.x()*b.x() + a.y()*b.y() + a.z()*b.z();
}

static inline bool hint_unlikely(bool expr)
{
#if (GCC_VERSION >= 302) || (__INTEL_COMPILER >= 800) || defined(__clang__)
  return __builtin_expect(expr != 0, 0);
#else
  return expr;
#endif
}

DEFUN_DLD(cbiot, args, nargout,
          "help string")
{
#ifndef NDEBUG 
  if ( hint_unlikely(args.length() != 4) )  {
    octave_stdout << "Expected 4 arguments: vihat = cbiot(pa, pb, pc, xh)\n";
    return octave_value_list();  
  }
  #endif 

  Matrix pa = args(0).matrix_value();
  dim_vector dm = pa.dims();
  intptr_t nvx = dm(0);

#ifndef NDEBUG
  if ( hint_unlikely(dm(1) != 3) )  {
    octave_stdout << "pa must have size [nvx,3]\n";
    return octave_value_list();  
  }
#endif 
  
  Matrix pb = args(1).matrix_value();
  dm = pb.dims();

#ifndef NDEBUG
  if ( hint_unlikely(dm(0) != nvx or dm(1) != 3) )  {
    octave_stdout << "pb must have size [nvx,3]\n";
    return octave_value_list();  
  }
#endif 

  Matrix pc = args(2).matrix_value();
  dm = pc.dims();
  intptr_t ncp = dm(0);

#ifndef NDEBUG  
  if ( hint_unlikely(dm(1) != 3) )  {
    octave_stdout << "pc must have size [ncp,3]\n";
    return octave_value_list();  
  }
#endif 

  Matrix xh = args(3).matrix_value();
  dm = xh.dims();

#ifndef NDEBUG
  if ( hint_unlikely(dm(0) != nvx or dm(1) != 3) )  {
    octave_stdout << "xh must have size [nvx,3]\n";
    return octave_value_list();  
  }
#endif 

  Matrix vihat(nvx, 3*ncp);
  double *pvihat = vihat.fortran_vec();
  std::fill(pvihat, pvihat + 3*nvx*ncp, 0.0);

  PListView la( pa.fortran_vec(), nvx );
  PListView lb( pb.fortran_vec(), nvx );
  PListView lx( xh.fortran_vec(), nvx );
  PListView lc( pc.fortran_vec(), ncp );
  const double i4pi = 1.0 / (4.0*M_PI);
  intptr_t npairs = ncp*nvx;

#pragma omp parallel if(npairs >= 4096)
  {
    double tmp[4*3];
    PListView ltmp(tmp, 4);
    PRef a(ltmp, 0);
    PRef b(ltmp, 1);
    PRef axx(ltmp, 2);
    PRef bxx(ltmp, 3);
#pragma omp for schedule(static,256)
    for (intptr_t i=0; i<ncp; ++i) {
      PRef cp( lc, i );
      PListView lv(pvihat + 3*i*nvx, nvx);

#pragma clang loop unroll(enable)
      for (intptr_t j=0; j<nvx; ++j) {
        PRef pta(la, j);
        PRef ptb(lb, j);
        PRef pxh(lx, j);

        a.sub(cp, pta);
        axx.cross(a, pxh);
        double adx = dot(a, pxh);
        double sqa = dot(a, a);

        b.sub(cp, ptb);
        bxx.cross(b, pxh);
        double bdx = dot(b, pxh);
        double sqb = dot(b, b);

        PRef pv(lv, j);
        axx.scale( 1.0 / (sqa - sqrt(sqa)*adx) );
        bxx.scale( 1.0 / (sqb - sqrt(sqb)*bdx) );
        pv.sub(axx, bxx);

        // self-induction term
        const double eps = std::numeric_limits<double>::epsilon();
        double t1 = sqrt(sqa*sqb) + dot(a, b);
        if (fabs(t1) > eps) {
          axx.cross(a,b);
          double t2 = 1.0 / t1;
          pv.x() += axx.x() * t2;
          pv.y() += axx.y() * t2;
          pv.z() += axx.z() * t2;
        }

        // divide by 4*pi
        pv.scale( i4pi );
      }
    }
  }  

  return octave_value(vihat);
}

