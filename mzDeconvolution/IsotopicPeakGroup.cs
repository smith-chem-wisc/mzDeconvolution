using System;
using IO.Thermo;
using System.Collections.Generic;
using Chemistry;
using System.Linq;
using Spectra;

namespace mzDeconvolution
{
    internal class IsotopicPeakGroup
    {
        public List<Peak> peakList = new List<Peak>();
        public int charge;
        public int Count
        {
            get
            {
                return peakList.Count;
            }
        }
        public double MostIntenseMass
        {
            get
            {
                return peakList.OrderByDescending(b => b.Y).First().X.ToMass(charge);
            }

        }
        public override string ToString()
        {
            return "C: " + charge;
        }

        public IsotopicPeakGroup(ThermoMzPeak peak)
        {
            peakList.Add(peak);
            this.charge = peak.Charge;
        }

        internal bool AttemptToAddNextIsotopicPeak(ThermoMzPeak peak)
        {
            if (peak.Charge == charge &&
                peak.MZ.ToMass(charge) < (peakList.Last().X.ToMass(charge) + 1) + Program.tol && peak.MZ.ToMass(charge) > (peakList.Last().X.ToMass(charge) + 1) - Program.tol)
            {
                peakList.Add(peak);
                return true;
            }
            return false;
        }
    }
}