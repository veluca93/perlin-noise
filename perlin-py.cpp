#include "perlin.hpp"
#include <functional>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace pybind11::literals;

template <typename T, typename U>
std::function<T(U &, py::args)>
args_to_doublev(T (U::*fpr)(const std::vector<double> &)) {
  return [fpr](U &f, py::args args) {
    return (f.*fpr)(args.cast<std::vector<double>>());
  };
};

template <typename T, typename U>
std::function<T(U &, py::args)>
args_to_rangev(T (U::*fpr)(const std::vector<PerlinNoiseFactory::Range> &)) {
  return [fpr](U &f, py::args args) {
    std::vector<PerlinNoiseFactory::Range> r;
    for (auto a : args) {
      auto pts = a.cast<std::vector<double>>();
      switch (pts.size()) {
      case 2:
        r.emplace_back(pts[0], pts[1]);
        break;
      case 3:
        r.emplace_back(pts[0], pts[1], pts[2], 1);
        break;
      case 4:
        r.emplace_back(pts[0], pts[1], pts[2], pts[3]);
        break;
      default:
        throw std::runtime_error("Invalid range specification, expected from 2 "
                                 "to 4 numbers but got " +
                                 std::to_string(pts.size()));
      }
    }
    return (f.*fpr)(r);
  };
};

PYBIND11_MODULE(cperlin, m) {
  py::class_<PerlinNoiseFactory>(m, "PerlinNoiseFactory")
      .def(py::init<>(&PerlinNoiseFactory::New), "dimension"_a, "octaves"_a = 1,
           "tile"_a = py::tuple{}, "unbias"_a = false, "seed"_a = 0)
      .def("get_plain_noise",
           args_to_doublev(&PerlinNoiseFactory::get_plain_noise))
      .def("__call__", args_to_doublev(&PerlinNoiseFactory::operator()))
      .def("get_byte", args_to_doublev(&PerlinNoiseFactory::get_byte))
      .def("get_range", args_to_rangev(&PerlinNoiseFactory::get_range))
      .def("get_range_as_bytes",
           args_to_rangev(&PerlinNoiseFactory::get_range_as_bytes));
}
