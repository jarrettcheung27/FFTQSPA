#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "BCJRQSPA.h"

namespace py = pybind11;

namespace {
int info_bits_len(const BCJRQSPA &codec) {
	return codec.get_block_length() * codec.m_degree;
}

int code_bits_len(const BCJRQSPA &codec) {
	return codec.get_total_length() * codec.m_degree;
}

void require_1d(const py::buffer_info &info, const char *name) {
	if (info.ndim != 1) {
		throw std::runtime_error(std::string(name) + " must be 1-D");
	}
}
}

PYBIND11_MODULE(fftqspa, m) {
	m.doc() = "FFTQSPA BiBo encoder/decoder bindings";

	py::class_<BCJRQSPA>(m, "BCJRQSPA")
		.def(py::init<std::string, int, std::string>(), py::arg("parity_filename"), py::arg("max_iteration"), py::arg("mapping_filename"))
		.def("info_bits_len", [](const BCJRQSPA &self) { return info_bits_len(self); })
		.def("code_bits_len", [](const BCJRQSPA &self) { return code_bits_len(self); })
		.def("encoder4bibo", [](BCJRQSPA &self, py::array_t<int, py::array::c_style | py::array::forcecast> info_bits) {
			auto in_info = info_bits.request();
			require_1d(in_info, "info_bits");

			const int expected = info_bits_len(self);
			if (in_info.size != expected) {
				throw std::runtime_error("info_bits length mismatch: expected " + std::to_string(expected));
			}

			py::array_t<int> code_bits(code_bits_len(self));
			auto out_info = code_bits.request();

			self.encoder4BiBo(static_cast<int *>(in_info.ptr), static_cast<int *>(out_info.ptr));
			return code_bits;
		}, py::arg("info_bits"), R"doc(
Encode binary info bits to binary codeword bits.

Parameters
----------
info_bits : 1-D int array
    Length must be info_bits_len(). Values should be 0/1.

Returns
-------
code_bits : 1-D int array
    Length is code_bits_len(). Values are 0/1.
)doc")
		.def("decode4bibo", [](BCJRQSPA &self, py::array_t<double, py::array::c_style | py::array::forcecast> rr_bits_prob) {
			auto in_info = rr_bits_prob.request();
			require_1d(in_info, "rr_bits_prob");

			const int expected = code_bits_len(self);
			if (in_info.size != expected) {
				throw std::runtime_error("rr_bits_prob length mismatch: expected " + std::to_string(expected));
			}

			py::array_t<int> decoded_bits(expected);
			auto out_info = decoded_bits.request();

			int iter = self.FFTQSPA4BiBo(static_cast<double *>(in_info.ptr), static_cast<int *>(out_info.ptr));
			return py::make_tuple(decoded_bits, iter);
		}, py::arg("rr_bits_prob"), R"doc(
Decode using FFTQSPA4BiBo.

Parameters
----------
rr_bits_prob : 1-D float array
    Length must be code_bits_len(). Each value is P(bit=0).

Returns
-------
(decoded_bits, iter) : (1-D int array, int)
    decoded_bits length is code_bits_len().
)doc")
		;
}
