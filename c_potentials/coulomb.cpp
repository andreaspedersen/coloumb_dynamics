/* File : coulomb.cpp */

#include "coulomb.h"

Coulomb::Coulomb()
{
    _nrAtoms = 0;
    setPrefactor(1.0);

    _xBox = 0.;
    _yBox = 0.;
    _zBox = 0.;
    _R0 = 0.;
    _energy = 0.;
    return;
}

Coulomb::~Coulomb()
{
    if (0 < _nrAtoms)
    {
        delete [] _positions;
        delete [] _forces;
        delete [] _charges;
    }
    return;
}

void Coulomb::setNrAtoms(int nrAtoms)
{
    if (0 < _nrAtoms)
    {
        delete [] _positions;
        delete [] _forces;
        delete [] _charges;
    }
    _positions = new double[3*nrAtoms];
    _forces = new double[3*nrAtoms];
    _charges = new double[nrAtoms];
    _nrAtoms = nrAtoms;
    return;
}

void Coulomb::setBox(double xBox, double yBox, double zBox)
{
    _needUpdate = 1;
    _needUpdateSingle = 1;

    _xBox = xBox;
    _yBox = yBox;
    _zBox = zBox;
    return;
}

void Coulomb::init(int nrAtoms)
{
    _needUpdate = 1;
    _needUpdateSingle = 1;

    setNrAtoms(nrAtoms);
    setPrefactor(1.0);
    return;
}

void Coulomb::setPrefactor(double prefactor)
{
    _needUpdate = 1;
    _needUpdateSingle = 1;

    _prefactor = prefactor;
    return;
}

void Coulomb::setR0(double R0)
{
    _needUpdate = 1;
    _needUpdateSingle = 1;

    _R0 = R0*R0;
    return;
}

void Coulomb::setPosition(int i, double x, double y, double z)
{
    _needUpdate = 1;
    _needUpdateSingle = 1;

    _positions[3*i  ] = x;
    _positions[3*i+1] = y;
    _positions[3*i+2] = z;

    _charges[i] = 1;

    return;
}

void Coulomb::setPositionAndCharge(int i, double x, double y, double z, double charge)
{
    _needUpdate = 1;
    _needUpdateSingle = 1;

    _positions[3*i  ] = x;
    _positions[3*i+1] = y;
    _positions[3*i+2] = z;

    _charges[i] = charge;

    return;
}

double Coulomb::getPosition(int i, int xyz)
{
    return _positions[3*i+xyz];
}

double Coulomb::getForce(int i, int xyz)
{
    if(_needUpdate)
    {
        _needUpdate = 0;
        _update();
    }
    return _forces[3*i+xyz];
}

double Coulomb::getForceSingle(int i, int xyz)
{
    if((_needUpdateSingle) || (_i != i))
    {
        _needUpdateSingle = 0;
        _updateSingle(i);
        _i = i;
    }
    return _forcesSingle[xyz];
}

double Coulomb::getEnergy()
{
    if(_needUpdate)
    {
        _needUpdate = 0;
        _update();
    }
    return _energy;
}

void Coulomb::_update()
{
    int index1 = 0;
    int index2 = 0;

    double distX = 0;
    double distY = 0;
    double distZ = 0;
    double normDistX = 0;
    double normDistY = 0;
    double normDistZ = 0;
    double distR = 0;
    double distR2 = 0;
    double f = 0;

    _energy = 0.;
    for(index1 = 0; index1 < _nrAtoms; index1++)
    {
        _forces[3*index1  ] = 0;
        _forces[3*index1+1] = 0;
        _forces[3*index1+2] = 0;
    }

    for(index1 = 0; index1 < _nrAtoms; index1++)
    {
        for(index2 = index1+1; index2 < _nrAtoms; index2++)
        {
            distX = _positions[3*index1  ] - _positions[3*index2  ];
            distY = _positions[3*index1+1] - _positions[3*index2+1];
            distZ = _positions[3*index1+2] - _positions[3*index2+2];

            if(_xBox != 0)
            {
                if(distX < -_xBox)
                    distX = distX + 2*_xBox;
                else if(distX > _xBox)
                    distX = distX - 2*_xBox;
            }
            if(_yBox != 0)
            {
                if(distY < -_yBox)
                    distY = distY + 2*_yBox;
                else if(distY > _yBox)
                    distY = distY - 2*_yBox;
            }
            if(_zBox != 0)
            {
                if(distZ < -_zBox)
                    distZ = distZ + 2*_zBox;
                else if(distZ > _zBox)
                    distZ = distZ - 2*_zBox;
            }
            distR2 = distX*distX + distY*distY + distZ*distZ;
            f = _prefactor / distR2 * _charges[index1] * _charges[index2];
            distR = sqrt(distR2);
            normDistX = distX / distR;
            normDistY = distY / distR;
            normDistZ = distZ / distR;

            // to prevent singularity at r = 0
            if(_R0 != 0.)
            {
                distR2 = distX*distX + distY*distY + distZ*distZ + _R0;
                f = _prefactor/distR2;
            }
            _forces[3*index1  ] = _forces[3*index1  ] + f*normDistX;
            _forces[3*index1+1] = _forces[3*index1+1] + f*normDistY;
            _forces[3*index1+2] = _forces[3*index1+2] + f*normDistZ;
            _forces[3*index2  ] = _forces[3*index2  ] - f*normDistX;
            _forces[3*index2+1] = _forces[3*index2+1] - f*normDistY;
            _forces[3*index2+2] = _forces[3*index2+2] - f*normDistZ;

            _energy = _energy + _prefactor / distR;
        }
    }
    return;
}
void Coulomb::_updateSingle(int i)
{
    int index1 = 0;
    double distX = 0;
    double distY = 0;
    double distZ = 0;
    double distR = 0;
    double distR2 = 0;
    double f = 0;

    _forcesSingle[0] = 0;
    _forcesSingle[1] = 0;
    _forcesSingle[2] = 0;

    for(index1 = 0; index1 < _nrAtoms; index1++)
    {
        if(index1 != i)
        {
            distX = _positions[3*index1  ] - _positions[3*i  ];
            distY = _positions[3*index1+1] - _positions[3*i+1];
            distZ = _positions[3*index1+2] - _positions[3*i+2];

            if(_xBox != 0)
            {
                if(distX < -_xBox)
                    distX = distX + 2*_xBox;
                else if(distX > _xBox)
                    distX = distX - 2*_xBox;
            }
            if(_yBox != 0)
            {
                if(distY < -_yBox)
                    distY = distY + 2*_yBox;
                else if(distY > _yBox)
                    distY = distY - 2*_yBox;
            }
            if(_zBox != 0)
            {
                if(distZ < -_zBox)
                    distZ = distZ + 2*_zBox;
                else if(distZ > _zBox)
                    distZ = distZ - 2*_zBox;
            }
            distR2 = distX*distX + distY*distY + distZ*distZ;
            f = _prefactor / distR2 * _charges[index1] * _charges[i];
            distR = sqrt(distR2);
            distX = distX / distR;
            distY = distY / distR;
            distZ = distZ / distR;

            // to prevent singularity at r = 0
            if(_R0 != 0.)
            {
                distR2 = distX*distX + distY*distY + distZ*distZ + _R0;
                f = _prefactor/distR2;
            }
            _forcesSingle[0] = _forcesSingle[0] - f*distX;
            _forcesSingle[1] = _forcesSingle[1] - f*distY;
            _forcesSingle[2] = _forcesSingle[2] - f*distZ;
        }
    }
    return;
}