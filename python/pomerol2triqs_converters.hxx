// DO NOT EDIT
// Generated automatically using libclang using the command :
// c++2py.py ../c++/pomerol_ed.hpp -I../../../local/pomerol/include -I/usr/include/eigen3 -I../c++ -p -mpytriqs.applications.impurity_solvers.pomerol2triqs -o pomerol2triqs --moduledoc "TRIQS wrapper around Pomerol ED library"


// --- C++ Python converter for g2_iw_freq_vec_params_t
#include <triqs/python_tools/converters/vector.hpp>
#include <triqs/python_tools/converters/string.hpp>
#include <algorithm>

namespace triqs { namespace py_tools {

template <> struct py_converter<g2_iw_freq_vec_params_t> {
 static PyObject *c2py(g2_iw_freq_vec_params_t const & x) {
  PyObject * d = PyDict_New();
  PyDict_SetItemString( d, "gf_struct"       , convert_to_python(x.gf_struct));
  PyDict_SetItemString( d, "beta"            , convert_to_python(x.beta));
  PyDict_SetItemString( d, "channel"         , convert_to_python(x.channel));
  PyDict_SetItemString( d, "three_freqs"     , convert_to_python(x.three_freqs));
  PyDict_SetItemString( d, "vec_four_indices", convert_to_python(x.vec_four_indices));
  return d;
 }

 template <typename T, typename U> static void _get_optional(PyObject *dic, const char *name, T &r, U const &init_default) {
  if (PyDict_Contains(dic, pyref::string(name)))
   r = convert_from_python<T>(PyDict_GetItemString(dic, name));
  else
   r = init_default;
 }

 template <typename T> static void _get_optional(PyObject *dic, const char *name, T &r) {
  if (PyDict_Contains(dic, pyref::string(name)))
   r = convert_from_python<T>(PyDict_GetItemString(dic, name));
  else
   r = T{};
 }

 static g2_iw_freq_vec_params_t py2c(PyObject *dic) {
  g2_iw_freq_vec_params_t res;
  res.gf_struct = convert_from_python<gf_struct_t>(PyDict_GetItemString(dic, "gf_struct"));
  res.beta = convert_from_python<double>(PyDict_GetItemString(dic, "beta"));
  _get_optional(dic, "channel"         , res.channel            ,PH);
  res.three_freqs = convert_from_python<three_freqs_t>(PyDict_GetItemString(dic, "three_freqs"));
  res.vec_four_indices = convert_from_python<std::vector<four_indices_t>>(PyDict_GetItemString(dic, "vec_four_indices"));
  return res;
 }

 template <typename T>
 static void _check(PyObject *dic, std::stringstream &fs, int &err, const char *name, const char *tname) {
  if (!convertible_from_python<T>(PyDict_GetItemString(dic, name), false))
   fs << "\n" << ++err << " The parameter " << name << " does not have the right type : expecting " << tname
      << " in C++, but got '" << PyDict_GetItemString(dic, name)->ob_type->tp_name << "' in Python.";
 }

 template <typename T>
 static void _check_mandatory(PyObject *dic, std::stringstream &fs, int &err, const char *name, const char *tname) {
  if (!PyDict_Contains(dic, pyref::string(name)))
   fs << "\n" << ++err << " Mandatory parameter " << name << " is missing.";
  else _check<T>(dic,fs,err,name,tname);
 }

 template <typename T>
 static void _check_optional(PyObject *dic, std::stringstream &fs, int &err, const char *name, const char *tname) {
  if (PyDict_Contains(dic, pyref::string(name))) _check<T>(dic, fs, err, name, tname);
 }

 static bool is_convertible(PyObject *dic, bool raise_exception) {
  if (dic == nullptr or !PyDict_Check(dic)) {
   if (raise_exception) { PyErr_SetString(PyExc_TypeError, "The function must be called with named arguments");}
   return false;
  }
  std::stringstream fs, fs2; int err=0;

#ifndef TRIQS_ALLOW_UNUSED_PARAMETERS
  std::vector<std::string> ks, all_keys = {"gf_struct","beta","channel","three_freqs","vec_four_indices"};
  pyref keys = PyDict_Keys(dic);
  if (!convertible_from_python<std::vector<std::string>>(keys, true)) {
   fs << "\nThe dict keys are not strings";
   goto _error;
  }
  ks = convert_from_python<std::vector<std::string>>(keys);
  for (auto & k : ks)
   if (std::find(all_keys.begin(), all_keys.end(), k) == all_keys.end())
    fs << "\n"<< ++err << " The parameter '" << k << "' is not recognized.";
#endif

  _check_mandatory<gf_struct_t                >(dic, fs, err, "gf_struct"       , "gf_struct_t");
  _check_mandatory<double                     >(dic, fs, err, "beta"            , "double");
  _check_optional <pomerol2triqs::channel_t   >(dic, fs, err, "channel"         , "pomerol2triqs::channel_t");
  _check_mandatory<three_freqs_t              >(dic, fs, err, "three_freqs"     , "three_freqs_t");
  _check_mandatory<std::vector<four_indices_t>>(dic, fs, err, "vec_four_indices", "std::vector<four_indices_t>");
  if (err) goto _error;
  return true;

 _error:
   fs2 << "\n---- There " << (err > 1 ? "are " : "is ") << err<< " error"<<(err >1 ?"s" : "")<< " in Python -> C++ transcription for the class g2_iw_freq_vec_params_t\n" <<fs.str();
   if (raise_exception) PyErr_SetString(PyExc_TypeError, fs2.str().c_str());
  return false;
 }
};

}}


// --- C++ Python converter for g2_iw_freq_box_params_t
#include <triqs/python_tools/converters/vector.hpp>
#include <triqs/python_tools/converters/string.hpp>
#include <algorithm>

namespace triqs { namespace py_tools {

template <> struct py_converter<g2_iw_freq_box_params_t> {
 static PyObject *c2py(g2_iw_freq_box_params_t const & x) {
  PyObject * d = PyDict_New();
  PyDict_SetItemString( d, "gf_struct", convert_to_python(x.gf_struct));
  PyDict_SetItemString( d, "beta"     , convert_to_python(x.beta));
  PyDict_SetItemString( d, "channel"  , convert_to_python(x.channel));
  PyDict_SetItemString( d, "n_b"      , convert_to_python(x.n_b));
  PyDict_SetItemString( d, "n_f"      , convert_to_python(x.n_f));
  PyDict_SetItemString( d, "index1"   , convert_to_python(x.index1));
  PyDict_SetItemString( d, "index2"   , convert_to_python(x.index2));
  PyDict_SetItemString( d, "index3"   , convert_to_python(x.index3));
  PyDict_SetItemString( d, "index4"   , convert_to_python(x.index4));
  return d;
 }

 template <typename T, typename U> static void _get_optional(PyObject *dic, const char *name, T &r, U const &init_default) {
  if (PyDict_Contains(dic, pyref::string(name)))
   r = convert_from_python<T>(PyDict_GetItemString(dic, name));
  else
   r = init_default;
 }

 template <typename T> static void _get_optional(PyObject *dic, const char *name, T &r) {
  if (PyDict_Contains(dic, pyref::string(name)))
   r = convert_from_python<T>(PyDict_GetItemString(dic, name));
  else
   r = T{};
 }

 static g2_iw_freq_box_params_t py2c(PyObject *dic) {
  g2_iw_freq_box_params_t res;
  res.gf_struct = convert_from_python<gf_struct_t>(PyDict_GetItemString(dic, "gf_struct"));
  res.beta = convert_from_python<double>(PyDict_GetItemString(dic, "beta"));
  _get_optional(dic, "channel"  , res.channel     ,PH);
  res.n_b = convert_from_python<int>(PyDict_GetItemString(dic, "n_b"));
  res.n_f = convert_from_python<int>(PyDict_GetItemString(dic, "n_f"));
  res.index1 = convert_from_python<indices_t>(PyDict_GetItemString(dic, "index1"));
  res.index2 = convert_from_python<indices_t>(PyDict_GetItemString(dic, "index2"));
  res.index3 = convert_from_python<indices_t>(PyDict_GetItemString(dic, "index3"));
  res.index4 = convert_from_python<indices_t>(PyDict_GetItemString(dic, "index4"));
  return res;
 }

 template <typename T>
 static void _check(PyObject *dic, std::stringstream &fs, int &err, const char *name, const char *tname) {
  if (!convertible_from_python<T>(PyDict_GetItemString(dic, name), false))
   fs << "\n" << ++err << " The parameter " << name << " does not have the right type : expecting " << tname
      << " in C++, but got '" << PyDict_GetItemString(dic, name)->ob_type->tp_name << "' in Python.";
 }

 template <typename T>
 static void _check_mandatory(PyObject *dic, std::stringstream &fs, int &err, const char *name, const char *tname) {
  if (!PyDict_Contains(dic, pyref::string(name)))
   fs << "\n" << ++err << " Mandatory parameter " << name << " is missing.";
  else _check<T>(dic,fs,err,name,tname);
 }

 template <typename T>
 static void _check_optional(PyObject *dic, std::stringstream &fs, int &err, const char *name, const char *tname) {
  if (PyDict_Contains(dic, pyref::string(name))) _check<T>(dic, fs, err, name, tname);
 }

 static bool is_convertible(PyObject *dic, bool raise_exception) {
  if (dic == nullptr or !PyDict_Check(dic)) {
   if (raise_exception) { PyErr_SetString(PyExc_TypeError, "The function must be called with named arguments");}
   return false;
  }
  std::stringstream fs, fs2; int err=0;

#ifndef TRIQS_ALLOW_UNUSED_PARAMETERS
  std::vector<std::string> ks, all_keys = {"gf_struct","beta","channel","n_b","n_f","index1","index2","index3","index4"};
  pyref keys = PyDict_Keys(dic);
  if (!convertible_from_python<std::vector<std::string>>(keys, true)) {
   fs << "\nThe dict keys are not strings";
   goto _error;
  }
  ks = convert_from_python<std::vector<std::string>>(keys);
  for (auto & k : ks)
   if (std::find(all_keys.begin(), all_keys.end(), k) == all_keys.end())
    fs << "\n"<< ++err << " The parameter '" << k << "' is not recognized.";
#endif

  _check_mandatory<gf_struct_t             >(dic, fs, err, "gf_struct", "gf_struct_t");
  _check_mandatory<double                  >(dic, fs, err, "beta"     , "double");
  _check_optional <pomerol2triqs::channel_t>(dic, fs, err, "channel"  , "pomerol2triqs::channel_t");
  _check_mandatory<int                     >(dic, fs, err, "n_b"      , "int");
  _check_mandatory<int                     >(dic, fs, err, "n_f"      , "int");
  _check_mandatory<indices_t               >(dic, fs, err, "index1"   , "indices_t");
  _check_mandatory<indices_t               >(dic, fs, err, "index2"   , "indices_t");
  _check_mandatory<indices_t               >(dic, fs, err, "index3"   , "indices_t");
  _check_mandatory<indices_t               >(dic, fs, err, "index4"   , "indices_t");
  if (err) goto _error;
  return true;

 _error:
   fs2 << "\n---- There " << (err > 1 ? "are " : "is ") << err<< " error"<<(err >1 ?"s" : "")<< " in Python -> C++ transcription for the class g2_iw_freq_box_params_t\n" <<fs.str();
   if (raise_exception) PyErr_SetString(PyExc_TypeError, fs2.str().c_str());
  return false;
 }
};

}}