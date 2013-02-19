
#ifndef COULOMB_H
#define COULOMB_H
/* File : coulomb.h */

#include <math.h>

class Coulomb{
public:
    Coulomb();
    ~Coulomb();
    void init(int nrAtoms);
    void setPrefactor(double prefactor);
	void setR0(double R0);
    void setNrAtoms(int nrAtoms);
    void setBox(double x, double y, double z);
    void setPosition(int i, double x, double y, double z);
    void setPositionAndCharge(int i, double x, double y, double z, double charge);
    double getPosition(int i, int xyz);
    double getForce(int i, int xyz);
    double getForceSingle(int i, int xyz);
	double getEnergy();

private:
    void _update();
    void _updateSingle(int i);

    double *_positions;
    double *_charges;

    int _nrAtoms;
    double _xBox;
    double _yBox;
    double _zBox;

	double _energy;

	int _needUpdate;
    double *_forces;

	int _i;
	int _needUpdateSingle;
    double _forcesSingle[3];

    double _prefactor;
	double _R0;
};

#endif