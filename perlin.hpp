#ifndef PERLIN_HPP
#define PERLIN_HPP
#include <cassert>
#include <functional>
#include <iostream>
#include <math.h>
#include <memory>
#include <random>
#include <stdexcept>
#include <unordered_map>
#include <vector>

template <typename T> class HashContainer {
public:
  std::size_t operator()(const T &vec) const {
    std::size_t ret = 0;
    for (auto i : vec) {
      ret ^= std::hash<decltype(i)>()(i) + 0x9e3779b9 + (ret << 6) + (ret >> 2);
    }
    return ret;
  }
};

class PerlinNoiseFactory {
public:
  struct Range {
    double scale;
    size_t start;
    size_t stop;
    size_t step;
    size_t sz() const { return (stop - start) / step; }
    Range(double scale, size_t v) : Range(scale, v, v + 1, 1) {}
    Range(double scale, size_t start, size_t stop, size_t step = 1)
        : scale(scale), start(start), stop(stop), step(step) {}
    Range() : Range(0.0, 0, 0, 0) {}
  };

  virtual double get_plain_noise(const std::vector<double> &point) = 0;
  virtual double operator()(const std::vector<double> &point) = 0;
  virtual uint8_t get_byte(const std::vector<double> &point) = 0;
  virtual std::vector<double> get_range(const std::vector<Range> &range) = 0;
  virtual std::vector<uint8_t>
  get_range_as_bytes(const std::vector<Range> &range) = 0;

  static std::unique_ptr<PerlinNoiseFactory>
  New(int dimension, int octaves = 1, const std::vector<size_t> &tile = {},
      bool unbias = false, uint32_t seed = 0);

protected:
  int dimension_;
  int octaves_;
  std::vector<size_t> tile_;
  bool unbias_;
  double scale_factor_;
  std::mt19937_64 rng_;
};

template <int dim_, typename T> struct storage_for_ {
  using type = std::array<T, dim_>;
};
template <typename T> struct storage_for_<0, T> {
  using type = std::vector<T>;
};

template <typename storage_t> storage_t cast_(const std::vector<double> &v) {
  storage_t ret = {};
  for (size_t i = 0; i < v.size(); i++)
    ret[i] = v[i];
  return ret;
}

template <>
std::vector<double> cast_<std::vector<double>>(const std::vector<double> &v) {
  return v;
}

template <typename T, size_t x>
void resize_(std::array<T, x> &v, size_t size) {}

template <typename T> void resize_(std::vector<T> &v, size_t size) {
  v.resize(size);
}

template <typename T, typename U> class Map {
public:
  Map(std::function<U()> generate) : generate_(generate) {}
  U &Get(const T &key) {
    EnsurePresent(key);
    return mp_.at(key);
  }
  void EnsurePresent(const T &key) {
    if (mp_.count(key))
      return;
    mp_.emplace(key, generate_());
  }

private:
  std::unordered_map<T, U, HashContainer<T>> mp_;
  std::function<U()> generate_;
};

template <typename T, size_t sz, typename U> class Map<std::array<T, sz>, U> {
public:
  Map(std::function<U()> generate) : generate_(generate) {}
  U &Get(const std::array<T, sz> &key) { return Get(key, sz - 1); }
  void EnsurePresent(const std::array<T, sz> &key) {
    return EnsurePresent(key, sz - 1);
  }
  template <size_t sz2> U &Get(const std::array<T, sz2> &key, size_t pos) {
    static_assert(sz2 >= sz, "Invalid call to Get");
    EnsurePresent(key, pos);
    return mp_[key[pos]].Get(key, pos - 1);
  }
  template <size_t sz2>
  void EnsurePresent(const std::array<T, sz2> &key, size_t pos) {
    static_assert(sz2 >= sz, "Invalid call to EnsurePresent");
    while (mp_.size() <= key[pos]) {
      mp_.emplace_back(generate_);
    }
    return mp_[key[pos]].EnsurePresent(key, pos - 1);
  }

private:
  std::vector<Map<std::array<T, sz - 1>, U>> mp_;
  std::function<U()> generate_;
};

template <typename T, typename U> class Map<std::array<T, 1>, U> {
public:
  Map(std::function<U()> generate) : generate_(generate) {}
  U &Get(const std::array<T, 1> &key) { return Get(key, 0); }
  void EnsurePresent(const std::array<T, 0> &key) {
    return EnsurePresent(key, 0);
  }
  template <size_t sz2> U &Get(const std::array<T, sz2> &key, size_t pos) {
    static_assert(sz2 >= 1, "Invalid call to Get");
    EnsurePresent(key, pos);
    return mp_[key[pos]];
  }
  template <size_t sz2>
  void EnsurePresent(const std::array<T, sz2> &key, size_t pos) {
    static_assert(sz2 >= 1, "Invalid call to EnsurePresent");
    while (mp_.size() <= key[pos]) {
      mp_.emplace_back(generate_());
    }
  }

private:
  std::vector<U> mp_;
  std::function<U()> generate_;
};

template <int dim_> class PerlinNoiseFactoryImpl : public PerlinNoiseFactory {
  using point_t = typename storage_for_<dim_, double>::type;
  using grid_t = typename storage_for_<dim_, size_t>::type;

public:
  PerlinNoiseFactoryImpl(int dimension, int octaves,
                         const std::vector<size_t> &tile, bool unbias,
                         uint32_t seed)
      : dimension_(dimension), octaves_(octaves), tile_(tile), unbias_(unbias),
        rng_(seed) {
    if (dim_ != 0 && dimension != dim_) {
      throw std::runtime_error("Invalid dimension for implementation");
    }
    if (dimension < 1)
      throw std::runtime_error("Invalid dimension number");
    if (tile.size() == 0) {
      for (int _ = 0; _ < dimension; _++)
        tile_.push_back(0);
    }
    if (tile.size() != (size_t)dimension)
      throw std::runtime_error("invalid tile value");
    scale_factor_ = 2 * pow(dimension, -.5);
  }

  double get_plain_noise(const std::vector<double> &point) override {
    if (point.size() != (size_t)dimension()) {
      throw std::runtime_error("Expected " + std::to_string(dimension()) +
                               " values, got " + std::to_string(point.size()));
    }
    return get_plain_noise_inner(cast_<point_t>(point));
  }

  double operator()(const std::vector<double> &point) override {
    return call(cast_<point_t>(point));
  }

  uint8_t get_byte(const std::vector<double> &point) override {
    return double_to_byte((*this)(point));
  }

  std::vector<double> get_range(const std::vector<Range> &range) override {
    if (range.size() != (size_t)dimension()) {
      throw std::runtime_error("Expected " + std::to_string(dimension()) +
                               " values, got " + std::to_string(range.size()));
    }
    size_t sz = 1;
    for (auto r : range) {
      sz *= r.sz();
    }
    std::vector<double> res(sz);
    if (sz == 0)
      return res;
    std::vector<size_t> point(dimension());
    for (int i = 0; i < dimension(); i++) {
      point[i] = range[i].start;
    }
    point_t scaled_point;
    resize_(scaled_point, dimension());
    int dim = 0;
    size_t cur = 0;
    while (dim < dimension()) {
      for (int i = 0; i < dimension(); i++) {
        scaled_point[i] = point[i] / range[i].scale;
      }
      res[cur] = call(scaled_point);
      cur++;
      while (dim < dimension()) {
        point[dim] += range[dim].step;
        if (point[dim] >= range[dim].stop) {
          point[dim] = range[dim].start;
          dim++;
        } else {
          dim = 0;
          break;
        }
      }
    }
    assert(cur == sz);
    return res;
  }

  std::vector<uint8_t>
  get_range_as_bytes(const std::vector<Range> &range) override {
    auto temp = get_range(range);
    std::vector<uint8_t> res(temp.size());
    for (size_t i = 0; i < temp.size(); i++) {
      res[i] = double_to_byte(temp[i]);
    }
    return res;
  }

private:
  int dimension_;

  int dimension() const {
    if (dim_)
      return dim_;
    return dimension_;
  }

  int octaves_;
  std::vector<size_t> tile_;
  bool unbias_;
  double scale_factor_;
  std::mt19937_64 rng_;
  Map<grid_t, point_t> gradient_{[this]() { return generate_gradient(); }};
  static uint8_t double_to_byte(double v) { return (v + 1) / 2 * 255 + 0.5; }

  static double smoothstep(double t) { return t * t * (3. - 2. * t); }

  static double lerp(double t, double a, double b) { return a + t * (b - a); }

  point_t generate_gradient() {
    if (dimension() == 1) {
      return {std::uniform_real_distribution<>(-1, 1)(rng_)};
    }
    std::normal_distribution<> dist(0, 1);
    point_t ret;
    resize_(ret, dimension());
    double norm = 0.0;
    for (int i = 0; i < dimension(); i++) {
      ret[i] = dist(rng_);
      norm += ret[i] * ret[i];
    }
    norm = pow(norm, 0.5);
    for (int i = 0; i < dimension(); i++) {
      ret[i] /= norm;
    }
    return ret;
  }

  double get_plain_noise_inner(point_t point, int exp = 1) {
    grid_t base_point;
    resize_(base_point, dimension());
    double result = 0.0;
    for (int i = 0; i < dimension(); i++) {
      point[i] *= 1ULL << exp;
      base_point[i] = point[i];
      if (tile_[i]) {
        base_point[i] %= tile_[i] << exp;
      }
    }
    static auto interpolate_helper = point;
    for (uint64_t set = 0; set < (1ULL << (uint64_t)dimension()); set++) {
      for (int i = 0; i < dimension(); i++) {
        if (set & (1 << i))
          base_point[i] += 1;
      }
      auto mod_base_point = base_point;
      for (int i = 0; i < dimension(); i++) {
        if (tile_[i] && mod_base_point[i] == tile_[i] << exp) {
          mod_base_point[i] = 0;
        }
      }
      const auto &gradient = gradient_.Get(mod_base_point);
      double dot = 0.0;
      for (int i = 0; i < dimension(); i++) {
        dot += gradient[i] * (point[i] - base_point[i]);
      }
      // do interpolation
      for (int i = 0; i < dimension(); i++) {
        if (set & (1 << i)) {
          double s = smoothstep(point[i] + 1 - base_point[i]);
          dot = lerp(s, interpolate_helper[i], dot);
        } else {
          // Save the result to interpolate it later.
          interpolate_helper[i] = dot;
          break;
        }
      }
      if (set == (1ULL << (uint64_t)dimension()) - 1) {
        // At the last iteration, dot is the result of all
        // the interpolations (as all bits in set are 1).
        result = dot;
      }
      for (int i = 0; i < dimension(); i++) {
        if (set & (1 << i))
          base_point[i] -= 1;
      }
    }
    return result * scale_factor_;
  }
  double call(const point_t &point) {
    double ret = 0;
    for (int o = 0; o < octaves_; o++) {
      ret += get_plain_noise_inner(point, o) / (1ULL << o);
    }
    ret /= 2 - pow(2, 1 - octaves_);
    if (unbias_) {
      double r = (ret + 1) / 2;
      for (int _ = 0; _ < (octaves_ + 1) / 2; _++) {
        r = smoothstep(r);
      }
      ret = r * 2 - 1;
    }
    return ret;
  }
};

std::unique_ptr<PerlinNoiseFactory>
PerlinNoiseFactory::New(int dimension, int octaves,
                        const std::vector<size_t> &tile, bool unbias,
                        uint32_t seed) {
  if (dimension == 1)
    return std::make_unique<PerlinNoiseFactoryImpl<1>>(dimension, octaves, tile,
                                                       unbias, seed);
  if (dimension == 2)
    return std::make_unique<PerlinNoiseFactoryImpl<2>>(dimension, octaves, tile,
                                                       unbias, seed);
  if (dimension == 3)
    return std::make_unique<PerlinNoiseFactoryImpl<3>>(dimension, octaves, tile,
                                                       unbias, seed);
  if (dimension == 4)
    return std::make_unique<PerlinNoiseFactoryImpl<4>>(dimension, octaves, tile,
                                                       unbias, seed);
  if (dimension == 5)
    return std::make_unique<PerlinNoiseFactoryImpl<5>>(dimension, octaves, tile,
                                                       unbias, seed);
  return std::make_unique<PerlinNoiseFactoryImpl<0>>(dimension, octaves, tile,
                                                     unbias, seed);
}

#endif
