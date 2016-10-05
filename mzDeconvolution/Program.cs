using IO.Thermo;
using System;
using System.Linq;
using MassSpectrometry;
using System.Collections.Generic;

namespace mzDeconvolution
{
    class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine(args[0]);

            ThermoRawFile a = new ThermoRawFile(args[0]);
            a.Open();

            foreach (var ms1scan in a.Where(b => b.MsnOrder == 1))
            {
                Console.WriteLine("In MS1 scan " + ms1scan.ScanNumber);

                IEnumerable<CandidatePeakCollection> candidatePeakCollections = GetCandidatePeakCollections(ms1scan);

                foreach (var CandidatePeakCollection in candidatePeakCollections)
                {
                    Console.WriteLine(" Mass: " + GetMonoisotopicMass(CandidatePeakCollection));
                }
            }

            Console.Read();
        }

        private static double GetMonoisotopicMass(CandidatePeakCollection candidatePeakCollection)
        {
            return candidatePeakCollection.FirstMass();
        }

        private static IEnumerable<CandidatePeakCollection> GetCandidatePeakCollections(IMsDataScan<ThermoSpectrum> ms1scan)
        {
            List<CandidatePeakCollection> candidatePeakCollections = new List<CandidatePeakCollection>();
            Dictionary<int, List<ThermoMzPeak>> peaksByCharge = new Dictionary<int, List<ThermoMzPeak>>();
            foreach (var peak in ms1scan.MassSpectrum)
            {
                if (peak.Charge > 0)
                {
                    if (!peaksByCharge.ContainsKey(peak.Charge))
                        peaksByCharge.Add(peak.Charge, new List<ThermoMzPeak>());
                    peaksByCharge[peak.Charge].Add(peak);
                }
            }

            foreach (var charge_kvp in peaksByCharge)
            {
                List<CandidatePeakCollection> candidatePeakCollectionsForCharge = new List<CandidatePeakCollection>();
                int charge = charge_kvp.Key;

                foreach (var peak in charge_kvp.Value)
                {
                    bool peakAccepted = false;
                    foreach (var cpc in candidatePeakCollectionsForCharge)
                    {
                        if (cpc.AcceptPeak(peak))
                        {
                            peakAccepted = true;
                            break;
                        }
                    }
                    if (peakAccepted == false)
                        candidatePeakCollectionsForCharge.Add(new CandidatePeakCollection(peak, charge));
                }
                candidatePeakCollections.AddRange(candidatePeakCollectionsForCharge);
            }

            return candidatePeakCollections;
        }
    }
}
