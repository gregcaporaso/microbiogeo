#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "0.0.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

"""Module for biogeo statistical methods used in workflows."""

from numpy import isnan

class UnparsableLineError(Exception):
    def __init__(self, line):
        self.args = ("Encountered unparsable line: '%s'" % line,)


class UnparsableFileError(Exception):
    def __init__(self, method):
        self.args = ("Unable to parse %s results file." % method.DisplayName,)


class AbstractStatMethod(object):
    DirectoryName = None
    ResultsName = None
    DisplayName = None
    StatDisplayName = None

    def parse(self, results_f):
        raise NotImplementedError

    def parse_float(self, float_str, min_val=None, max_val=None,
                    suppress_nan_check=False):
        """Converts a float (as a string) into a float.

        Performs optional sanity checks to ensure the float is in
        [min_val, max_val].
        """
        try:
            result = float(float_str)
        except ValueError:
            raise ValueError("Could not convert float string '%s' to float."
                             % float_str)

        if not suppress_nan_check and isnan(result):
            raise TypeError("Encountered NaN instead of a valid float.")

        if (min_val is not None and result < min_val) or \
           (max_val is not None and result > max_val):
            raise ValueError(
                    "Float %.4f does not fall in valid range." % result)

        return result

    def __eq__(self, other):
        return type(self) == type(other)

    def __ne__(self, other):
        return not self.__eq__(other)


class QiimeStatMethod(AbstractStatMethod):
    def parse(self, results_f):
        for line in results_f:
            pass

        tokens = line.strip().split('\t')

        if len(tokens) != 4:
            raise UnparsableLineError(line)

        es, p_value = tokens[1:3]
        es = self.parse_float(es)

        if 'Too few iters to compute p-value' in p_value:
            raise UnparsableLineError(line)
        else:
            p_value = self.parse_float(p_value, 0, 1)

        return es, p_value


class Anosim(QiimeStatMethod):
    DirectoryName = 'anosim'
    ResultsName = 'anosim'
    DisplayName = 'ANOSIM'
    StatDisplayName = r'$R$'


class Permanova(QiimeStatMethod):
    DirectoryName = 'permanova'
    ResultsName = 'permanova'
    DisplayName = 'PERMANOVA'
    StatDisplayName = r'$F$'


class Adonis(AbstractStatMethod):
    DirectoryName = 'adonis'
    ResultsName = 'adonis'
    DisplayName = 'Adonis'
    StatDisplayName = r'$R^2$'

    def parse(self, results_f):
        for line in results_f:
            if line.startswith('qiime.data$map[[opts$category]]'):
                tokens = line.strip().split()

                # The format of the file changes if the result is significant
                # or not.
                if len(tokens) == 7:
                    es, p_value = tokens[-2:]
                elif len(tokens) == 8:
                    es, p_value = tokens[-3:-1]
                else:
                    raise UnparsableLineError(line)

                return (self.parse_float(es, 0, 1),
                        self.parse_float(p_value, 0, 1))

        raise UnparsableFileError(self)


class Mrpp(AbstractStatMethod):
    DirectoryName = 'mrpp'
    ResultsName = 'mrpp'
    DisplayName = 'MRPP'
    StatDisplayName = r'$A$'

    def parse(self, results_f):
        a_value = None
        p_value = None

        for line in results_f:
            if line.startswith('Chance corrected within-group agreement A:'):
                tokens = line.strip().split()

                if len(tokens) == 6:
                    a_value = self.parse_float(tokens[-1])
                else:
                    raise UnparsableLineError(line)
            elif line.startswith('Significance of delta:'):
                tokens = line.strip().split()

                if len(tokens) == 4:
                    p_value = self.parse_float(tokens[-1], 0, 1)
                else:
                    raise UnparsableLineError(line)

        if a_value is None or p_value is None:
            raise UnparsableFileError(self)

        return a_value, p_value


class Dbrda(AbstractStatMethod):
    DirectoryName = 'dbrda'
    ResultsName = 'dbrda'
    DisplayName = 'db-RDA'
    StatDisplayName = r'$R^2$'

    def parse(self, results_f):
        r2_value = None
        p_value = None

        for line in results_f:
            if line.startswith('Constrained'):
                tokens = line.strip().split()

                if len(tokens) == 4:
                    r2_value = self.parse_float(tokens[2], 0, 1)
                else:
                    raise UnparsableLineError(line)
            elif line.startswith('Significance:'):
                tokens = line.strip().split()

                if len(tokens) == 2:
                    p_value = self.parse_float(tokens[1], 0, 1)
                else:
                    raise UnparsableLineError(line)

        if r2_value is None or p_value is None:
            raise UnparsableFileError(self)

        return r2_value, p_value


class Permdisp(AbstractStatMethod):
    DirectoryName = 'permdisp'
    ResultsName = 'permdisp'
    DisplayName = 'PERMDISP'
    StatDisplayName = r'$F$'

    def parse(self, results_f):
        f_value = None
        p_value = None
        at_nonparametric_section = False

        for line in results_f:
            if line.startswith('No. of permutations:'):
                at_nonparametric_section = True
            elif line.startswith('Groups') and at_nonparametric_section:
                tokens = line.strip().split()

                if len(tokens) == 7 or len(tokens) == 8:
                    f_value = self.parse_float(tokens[4])
                    p_value = self.parse_float(tokens[6], 0, 1)
                else:
                    raise UnparsableLineError(line)

        if f_value is None or p_value is None:
            raise UnparsableFileError(self)

        return f_value, p_value


class Mantel(AbstractStatMethod):
    DirectoryName = 'mantel'
    ResultsName = 'mantel'
    DisplayName = 'Mantel'
    StatDisplayName = r'$r$'

    def parse(self, results_f):
        for line in results_f:
            pass

        tokens = line.strip().split('\t')

        if len(tokens) != 7:
            raise UnparsableLineError(line)

        es, p_value = tokens[3:5]
        es = self.parse_float(es, -1, 1, suppress_nan_check=True)
        p_value = self.parse_float(p_value, 0, 1)

        if isnan(es):
            es = 0.0
            p_value = 1.0

        return es, p_value


class PartialMantel(AbstractStatMethod):
    DirectoryName = 'partial_mantel'
    ResultsName = 'partial_mantel'
    DisplayName = 'Partial Mantel'
    StatDisplayName = r'$r$'

    def parse(self, results_f):
        for line in results_f:
            pass

        tokens = line.strip().split('\t')

        if len(tokens) != 8:
            raise UnparsableLineError(line)

        es, p_value = tokens[4:6]
        return self.parse_float(es, -1, 1), self.parse_float(p_value, 0, 1)


class MantelCorrelogram(AbstractStatMethod):
    DirectoryName = 'mantel_corr'
    DisplayName = 'Mantel Correlogram'


class MoransI(AbstractStatMethod):
    DirectoryName = 'morans_i'
    ResultsName = 'morans_i'
    DisplayName = 'Moran\'s I'
    StatDisplayName = r'$I$'

    def parse(self, results_f):
        es = None
        p_value = None
        es_next = False
        p_value_next = False

        for line in results_f:
            if line.startswith('$observed'):
                es_next = True
            elif line.startswith('$p.value'):
                p_value_next = True
            elif es_next:
                if len(line.strip().split()) != 2:
                    raise UnparsableLineError(line)

                es = self.parse_float(line.strip().split()[1], -1, 1)
                es_next = False
            elif p_value_next:
                if len(line.strip().split()) != 2:
                    raise UnparsableLineError(line)

                # Moran's I will sometimes calculate a p-value of 2.0. From
                # looking at the R code, it seems like this is a bug, and that
                # the p-value should be 1.0.
                p_value = self.parse_float(line.strip().split()[1], 0, 2)
                p_value = min(p_value, 1.0)
                p_value_next = False

        if es is None or p_value is None:
            raise UnparsableFileError(self)

        return es, p_value


class Best(AbstractStatMethod):
    DirectoryName = 'best'
    DisplayName = 'BEST'


class OrdinationCorrelation(AbstractStatMethod):
    ResultsName = 'ord_corr'

    def parse(self, results_f):
        for line in results_f:
            pass

        tokens = line.strip().split('\t')

        if len(tokens) != 3:
            raise UnparsableLineError(line)

        es = tokens[0]
        p_value = tokens[2]

        es = self.parse_float(es, -1, 1)

        if 'Too few iters to compute p-value' in p_value or p_value == 'N/A':
            raise UnparsableLineError(line)
        else:
            p_value = self.parse_float(p_value, 0, 1)

        return es, p_value


class PearsonOrdinationCorrelation(OrdinationCorrelation):
    DirectoryName = 'ord_corr_pearson'
    DisplayName = 'Ordination Correlation (Pearson)'
    StatDisplayName = r'$r$'


class SpearmanOrdinationCorrelation(OrdinationCorrelation):
    DirectoryName = 'ord_corr_spearman'
    DisplayName = 'Ordination Correlation (Spearman)'
    StatDisplayName = r'$\rho$'
