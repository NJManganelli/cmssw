#ifndef SimPPS_PPSPixelDigiProducer_ENERGY_DEPOSIT_UNIT_H
#define SimPPS_PPSPixelDigiProducer_ENERGY_DEPOSIT_UNIT_H

#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"

/**
 * Class which allows to "follow" an elementary charge in the silicon.
 * It basically defines a quantum of energy, with a position.
 */
class RPixEnergyDepositUnit {
public:
  RPixEnergyDepositUnit() : energy_(0), position_(0, 0, 0) {}

  RPixEnergyDepositUnit(double energy, double x, double y, double z) : energy_(energy), position_(x, y, z) {}

  RPixEnergyDepositUnit(double energy, const Local3DPoint& position) : energy_(energy), position_(position) {}

  inline double Energy() const { return energy_; }
  inline const Local3DPoint& Position() const { return position_; }

  inline void setEnergy(double e) { energy_ = e; }
  inline void setPosition(Local3DPoint p) { position_ = p; }

private:
  double energy_;
  Local3DPoint position_;
};

#endif  //SimPPS_PPSPixelDigiProducer_ENERGY_DEPOSIT_UNIT_H
