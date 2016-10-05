using System;
using IO.Thermo;
using System.Collections.Generic;
using Chemistry;
using System.Linq;

namespace mzDeconvolution
{
    internal class CandidatePeakCollection
    {
        private List<double> massList;
        private int charge;
        private static readonly double tol = 0.01;
        public CandidatePeakCollection(ThermoMzPeak peak, int charge)
        {
            massList = new List<double>();
            massList.Add(peak.MZ.ToMass(charge));
            this.charge = charge;
        }


        internal bool AcceptPeak(ThermoMzPeak peak)
        {
            if (peak.Charge == charge && peak.MZ.ToMass(charge) < (massList.Last() + 1) + tol && peak.MZ.ToMass(charge) > (massList.Last() + 1) - tol)
            {
                massList.Add(peak.MZ.ToMass(charge));
                return true;
            }
            return false;
        }

        internal double FirstMass()
        {
            return massList.First();
        }
    }
}