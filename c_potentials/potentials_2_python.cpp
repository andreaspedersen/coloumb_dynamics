/** @file
      Wrap Classes for Python.
      @author Andreas Pedersen, University of Reykjavik
      @date 2012
      */


#include "boost/python.hpp"

#include "coulomb.h"

using namespace std;
using namespace boost::python;

BOOST_PYTHON_MODULE(c_potentials)
{
	class_<Coulomb>("Coulomb_C")
    .def(init<>())
    .def("init", &Coulomb::init)
    .def("setBox", &Coulomb::setBox)
    .def("setNrAtoms", &Coulomb::setNrAtoms)
    .def("setPrefactor", &Coulomb::setPrefactor)
    .def("setR0", &Coulomb::setR0)
    .def("setPosition", &Coulomb::setPosition)
    .def("setPositionAndCharge", &Coulomb::setPositionAndCharge)
    .def("getForce", &Coulomb::getForce)
    .def("getForceSingle", &Coulomb::getForceSingle)
    .def("getEnergy", &Coulomb::getEnergy)
    ;
}
