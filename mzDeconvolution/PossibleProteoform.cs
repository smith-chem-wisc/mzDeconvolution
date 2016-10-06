﻿using Chemistry;
using System;
using System.Collections.Generic;
using System.Linq;

namespace mzDeconvolution
{
    internal class PossibleProteoform
    {
        public double EvidenceLevel { get; private set; }

        public HashSet<int> charges = new HashSet<int>();

        public HashSet<IsotopicPeakGroup> isotopicPeakGroups = new HashSet<IsotopicPeakGroup>();
        
        public double GetMonoisotopicMass()
        {
            double MostIntenseAverage = isotopicPeakGroups.Select(b => b.MostIntenseMass).Average();

            ChemicalFormula chemicalFormula = FindChemicalFormulaWithAverageMass(MostIntenseAverage);

            var ye = new IsotopicDistribution(chemicalFormula, 0.01, 0.001);
            double[] massesArray = ye.Masses.ToArray();
            double[] intensitiesArray = ye.Intensities.ToArray();
            Array.Sort(intensitiesArray, massesArray);
            double shiftFromMostIntenseToMono = massesArray.Last() - chemicalFormula.MonoisotopicMass;

            return MostIntenseMode() - shiftFromMostIntenseToMono;
        }

        private double MostIntenseMode()
        {
            Dictionary<double, List<double>> mostIntenses = new Dictionary<double, List<double>>();
            foreach (var ye in isotopicPeakGroups)
            {
                bool added = false;
                foreach (var kvp in mostIntenses)
                {
                    if (kvp.Key > ye.MostIntenseMass - 0.5 && kvp.Key < ye.MostIntenseMass + 0.5)
                    {
                        kvp.Value.Add(ye.MostIntenseMass);
                        added = true;
                        break;
                    }
                }
                if (!added)
                    mostIntenses.Add(ye.MostIntenseMass, new List<double>() { ye.MostIntenseMass });
            }
            int bestCount = 0;
            KeyValuePair<double, List<double>> bestKvp = new KeyValuePair<double, List<double>>();
            foreach (var kvp in mostIntenses)
            {
                if (kvp.Value.Count > bestCount)
                {
                    bestCount = kvp.Value.Count;
                    bestKvp = kvp;
                }
            }
            return bestKvp.Value.Average();
        }


        private ChemicalFormula FindChemicalFormulaWithAverageMass(double mostIntenseAverage)
        {
            const double averageC = 4.9384;
            const double averageH = 7.7583;
            const double averageO = 1.4773;
            const double averageN = 1.3577;
            const double averageS = 0.0417;

            ChemicalFormula chemicalFormula = new ChemicalFormula();

            double factor = (mostIntenseAverage / 98.123);

            chemicalFormula.Add("C", Convert.ToInt32(Math.Round(averageC * factor)));
            chemicalFormula.Add("H", Convert.ToInt32(Math.Round(averageH * factor)));
            chemicalFormula.Add("O", Convert.ToInt32(Math.Round(averageO * factor)));
            chemicalFormula.Add("N", Convert.ToInt32(Math.Round(averageN * factor)));
            chemicalFormula.Add("S", Convert.ToInt32(Math.Round(averageS * factor)));

            do
            {
                if (chemicalFormula.AverageMass > mostIntenseAverage + 33 && chemicalFormula.ContainsIsotopesOf(PeriodicTable.GetElement("S")))
                    chemicalFormula.Remove("S");
                if (chemicalFormula.AverageMass > mostIntenseAverage + 17 && chemicalFormula.ContainsIsotopesOf(PeriodicTable.GetElement("O")))
                    chemicalFormula.Remove("N");
                if (chemicalFormula.AverageMass > mostIntenseAverage + 15 && chemicalFormula.ContainsIsotopesOf(PeriodicTable.GetElement("N")))
                    chemicalFormula.Remove("O");
                if (chemicalFormula.AverageMass > mostIntenseAverage + 13 && chemicalFormula.ContainsIsotopesOf(PeriodicTable.GetElement("C")))
                    chemicalFormula.Remove("C");
                if (chemicalFormula.AverageMass > mostIntenseAverage + 1.001 && chemicalFormula.ContainsIsotopesOf(PeriodicTable.GetElement("H")))
                    chemicalFormula.Remove("H");

                if (chemicalFormula.AverageMass < mostIntenseAverage - 33)
                    chemicalFormula.Add("S");
                if (chemicalFormula.AverageMass < mostIntenseAverage - 17)
                    chemicalFormula.Add("N");
                if (chemicalFormula.AverageMass < mostIntenseAverage - 15)
                    chemicalFormula.Add("O");
                if (chemicalFormula.AverageMass < mostIntenseAverage - 13)
                    chemicalFormula.Add("C");
                if (chemicalFormula.AverageMass < mostIntenseAverage - 1.001)
                    chemicalFormula.Add("H");
            } while ((chemicalFormula.AverageMass > mostIntenseAverage + 1.001) || (chemicalFormula.AverageMass < mostIntenseAverage - 1.001));

            return chemicalFormula;
        }

        public PossibleProteoform(IsotopicPeakGroup isotopicPeakGroup1, IsotopicPeakGroup isotopicPeakGroup2, double score)
        {
            isotopicPeakGroups.Add(isotopicPeakGroup1);
            isotopicPeakGroups.Add(isotopicPeakGroup2);
            charges.Add(isotopicPeakGroup1.charge);
            charges.Add(isotopicPeakGroup2.charge);
            EvidenceLevel = score;
        }

        public override string ToString()
        {
            return "MonoisotopicMass: " + GetMonoisotopicMass() + " charges: " + string.Join(",", charges.OrderBy(b => b)) + " evidence level: " + EvidenceLevel;
        }

        internal void Add(IsotopicPeakGroup isotopicPeakGroup1, IsotopicPeakGroup isotopicPeakGroup2, double score)
        {
            isotopicPeakGroups.Add(isotopicPeakGroup1);
            isotopicPeakGroups.Add(isotopicPeakGroup2);
            charges.Add(isotopicPeakGroup1.charge);
            charges.Add(isotopicPeakGroup2.charge);
            EvidenceLevel += score;
        }
    }
}