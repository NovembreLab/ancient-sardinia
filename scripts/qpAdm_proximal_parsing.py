### Last modified 2019-06-17 by Chi-Chun

import os
import re
import itertools as it
import pandas as pd

# set working directory
print('current directory: ' + os.getcwd())
os.chdir('/project/jnovembre/jhmarcus/ancient-sardinia/output/qpAdm_rev/')
print('moved to directory: ' + os.getcwd()) 

### Parser Class
# log object
class qplog:
    '''
    Super class for all AdmixTools log files
    '''
                
    def print_anchor(self):
        print(self._anchor_line)
        
    def _parse_param_path(self):
        with open(self.logfile, 'r') as log:
            line = log.readline()
            info = line.strip().split(':')
            self.bin_path, _, self.par_path = [i.strip() for i in info]
    
    def _parse_version(self, line):
        with open(self.logfile, 'r') as log:
            for line in it.islice(log, line, line + 1):
                ver = (line.split(':')[1]).strip()
                self.version = int(ver)
    
    def __init__(self, logfile):
        self.logfile = logfile
        self.version = None
        self.data_info = {}
        self.parameter = {}
        self._parse_param_path()


class qpAdm_log(qplog):
    ''' 
    Class for Dstat test output
    '''
    def _get_anchor_line(self):
        self._anchor_line = {}
        with open(self.logfile, 'r') as log:
            for n, line in enumerate(log):
                if line.startswith('##PARAMETER NAME: VALUE'):
                    self._anchor_line['param_start_line'] = n + 1
#                 elif line.startswith('## qpAdm version:'):
#                     self._anchor_line['version_line'] = n
                elif line.startswith('left pops:'):
                    self._anchor_line['leftpop_start_line'] = n
#                 elif line.startswith('right pops:'):
#                     self._anchor_line['rightpop_start_line'] = n
#                 elif line.strip().startswith('0 '):
#                     self._anchor_line['sample_size_start_line'] = n
#                 elif line.startswith('jackknife block size:'):
#                      self._anchor_line['data_info_start_line'] = n
                elif line.startswith('## ncols:'):
                     self._anchor_line['coverage_start_line'] = n
                elif line.startswith('dof (jackknife):'):
                     self._anchor_line['used_data_info_start_line'] = n
                elif line.startswith('numsnps used:'):
                     self._anchor_line['used_data_info_end_line'] = n
                elif line.startswith('best coefficients:'):
                     self._anchor_line['summary_start_line'] = n
                elif line.startswith('coeffs:'):
                     self._anchor_line['summary_end_line'] = n
#                 elif line.startswith('## dscore::'):
#                     self._anchor_line['detail_start_line'] = n
                elif line.startswith('## end of run'):
                     self._anchor_line['end_line'] = n
    
    def _parse_param(self, start_line, end_line):
        with open(self.logfile, 'r') as log:
            for line in it.islice(log, start_line, end_line):
                param, val = line.strip().split(':')
                self.parameter[param] = val.strip()    
    
    def _parse_left_right_pops(self, start_line):
        pop_list = [] 
        with open(self.logfile, 'r') as log:
            for line in it.islice(log, start_line, None):
                if line.strip() != '':
                    pop_list.append(line.strip())
                else:
                    break
        return pop_list
        
    def _parse_sample_size(self, start_line, end_line):
        sample_sizes = []
        with open(self.logfile, 'r') as log:
            for line in it.islice(log, start_line, end_line):
                sizes = line.strip().split()
                sizes = [s.strip() for s in sizes]
                sample_sizes.append((sizes[1], sizes[2]))
        self.sample_size = pd.DataFrame(sample_sizes, columns = ['population', 'n'])
        
    def _parse_data_info(self, start_line, end_line):
        with open(self.logfile, 'r') as log:
            for line in it.islice(log, start_line, end_line):
                if line.startswith('snps'):
                    par1, val1, par2, val2 = line.strip().split()
                    self.data_info[par1] = int(val1.strip())
                    self.data_info[par2] = int(val2.strip())
                elif line.startswith('nrows'):
                    param, val = line.strip().split(':')
                    par1, par2 = param.split(',')
                    par2 = par2.strip()
                    val1, val2 = val.split()
                    self.data_info[(par1, par2)] = (int(val1), int(val2))
                else:
                    param, val = line.strip().split(':')
                    self.data_info[param] = float(val.strip())
    
    def _parse_coverage(self, start_line):
        coverage = []
        with open(self.logfile, 'r') as log:
            for line in it.islice(log, start_line, None):
                if not line.strip().startswith('coverage:'):
                    break
                line = re.split(r'\s+', line.strip())
                line[2] = int(line[2])
                coverage.append(line[1:3])
        self.coverage = pd.DataFrame(coverage, columns = ['population', 'coverage'])

    def get_f4_info(self, rank):
        '''rank, dof, chisq, tail, chisqdiff, taildiff
           A, B
        '''
        start_line = None
        with open(self.logfile, 'r') as log:
            # locate start line
            for n, line in enumerate(log):
                if line.startswith('f4rank: {}'.format(str(rank))):
                    start_line = n
                    print('start line {}'.format(start_line))
                    break
            assert start_line != None
        
        with open(self.logfile, 'r') as log:
            f4_info = {}
            # locate anchor lines and parse first line
            for i, line in enumerate(it.islice(log, start_line, None)):
                if i == 0:
                    line = line.strip().split()
                    param = line[::2]
                    val = line[1::2]
                    for j, p in enumerate(param):
                        f4_info[p] = val[j]
                elif line.strip() == 'B:':
                    i_B = i
                elif line.strip() == 'A:':
                    i_A = i
                elif line.strip() == '':
                    i_end = i
                    break
        
        with open(self.logfile, 'r') as log:
            # parse data
            B = []
            A = []
            for i, line in enumerate(it.islice(log, start_line, None)):

                temp_header = ['pop'] + [str(p) for p in range(rank)]
                #print(temp_header)
                if (i > i_B) & (i < i_A):
                    line = line.strip().split()
                    B.append(line)
                elif (i > i_A) & (i < i_end):
                    line = line.strip().split()
                    A.append(line)    
            f4_info['A'] = pd.DataFrame(A, columns = temp_header)
            f4_info['B'] = pd.DataFrame(B, columns = temp_header)
        
        return f4_info
    
    def _parse_summary(self, start_line, end_line):
        '''from (line) best coefficients after f4info to
           (line) coeffs before ##dscore:
        '''
        
        summary = {}
        table = []
        with open(self.logfile, 'r') as log:
            for i, line in enumerate(it.islice(log, start_line, end_line)):
                if line.strip().startswith(('best coefficients', 'Jackknife mean', 'std. errors')):
                    par = line.strip().split(':')[0]
                    line = line.strip().split(':')[1]
                    line = line.strip().split()
                    summary[par] = [float(l.strip()) for l in line]
                elif line.startswith('error covariance'):
                    part2_start_line = i + start_line
                elif line.startswith('summ:'):
                    part3_line = i + start_line
                elif line.strip().startswith('fixed pat'):
                    part4_start_line = i + start_line
                elif line.startswith('best pat:'):
                    p1 = line.strip().strip('best pat:').strip()
                    p1 = p1.split()[0:2]
                    p1[1] = float(p1[1])
                    if len(line.strip().split('chi(nested):')) > 1:
                        temp = line.strip().split('chi(nested):')[1].strip() 
                        if line.strip().endswith('infeasible'):
                            p2 = [float(temp.split()[0]), float((temp.split()[-2])), 'infeasible']
                        else:
                            p2 = [float(temp.split()[0]), float((temp.split()[-1])), '-']     
                    else:
                        p2 = [-1, -1, -1]
                    table.append(p1 + p2)
                elif line.startswith('coeffs:'):
                    pass
                else:
                    pass
                
        n_coeff = len(summary['std. errors'])
        
        if n_coeff > 1:
            summary['table'] = pd.DataFrame(table, columns = ['best pattern', 'p-val', 'chi(nested)', 'p-val for nested model', 'status'])
   
        pattern = []     
        cov = []
        
        with open(self.logfile, 'r') as log:
        # parse coefficient
            for i, line in enumerate(it.islice(log, part2_start_line + 1, None)):
                line = line.strip().split()
                cov.append([int(l) for l in line])           
                if i >= n_coeff - 1:
                    break
            summary['cov'] = cov
        
        with open(self.logfile, 'r') as log:
        # parse summ
            for i, line in enumerate(it.islice(log, part3_line, None)):
                
                if i == 0:
                    if not line.strip().endswith('...'):
                        line = line.strip().split()
                        line[2] = int(line[2])
                        line[3] = float(line[3])
                        line[4:(4 + n_coeff + 1)] = [float(l) for l in line[4:(4 + n_coeff + 1)]]
                        line[(4 + n_coeff + 1)::] = [int(l) for l in line[(4 + n_coeff + 1)::]]
                        summary['summ'] = line[1::]
                        break
                    else:
                        line = line.strip().split()[:-1]
                        line[2] = int(line[2])
                        line[4:(4 + n_coeff + 1)] = [float(l) for l in line[4:(4 + n_coeff + 1)]]
                        line[(4 + n_coeff + 1)::] = [int(l) for l in line[(4 + n_coeff + 1)::]]
                        line0 = line[1::]

                if i == 1:
                    line = line.strip().split()
                    line1 = [int(l) for l in line]
                    summary['summ'] = line0 + line1
                    
        if n_coeff > 1:
            with open(self.logfile, 'r') as log:
            # parse pattern table
                for i, line in enumerate(it.islice(log, part4_start_line, None)):
                    temp_header =  ['fixed pat', 'wt', 'dof', 'chisq', 'tail prob']
                    temp_header += ['coeff' + str(x) for x in range(1, n_coeff + 1)]
                    temp_header += ['status']

                    max_line = 2 ** n_coeff
                    if i == 0:
                        pass
                    elif (i < max_line):
                            line = line.strip().split()
                            line[1:3] = [int(l) for l in line[1:3]]
                            status = '-'
                            if line[-1] == 'infeasible':
                                line = line[:-1]
                                status = 'infeasible'
                            line[3:] = [float(l) for l in line[3:]]
                            line.append(status)
                            pattern.append(line)
                    summary['model_comparison'] = pd.DataFrame(pattern, columns = temp_header)
        
        self.summary = summary
            
    def _parse_details(self, start_line, end_line):
        '''
        dscore:: f_4(Base, Fit, Rbase, right2)
        gendstat:: f_4(Base, Fit, right1, right2)
        '''
        details = []
        dscores = []
        gendstats = []
        with open(self.logfile, 'r') as log:
            for line in it.islice(log, start_line, end_line):
                if line.startswith('details'):
                    line = line.strip().split()
                    line = [l.strip() for l in line]
                    details.append((line[1], line[2], float(line[3])))
                elif line.startswith('dscore'):
                    line = line.strip().split()
                    line = [l.strip() for l in line]
                    param = line[::2]
                    val = line[1::2]
                    dscores.append((val[0], float(val[1]), float(val[2])))
                elif line.startswith('gendstat:'):
                    line = line.strip().split()
                    line = [l.strip() for l in line]
                    gendstats.append((line[1], line[2], float(line[3])))
                else:
                    pass

            details = pd.DataFrame(details, columns = ['right1', 'right2', 'f4'])
            dscores = pd.DataFrame(dscores, columns = ['right1', 'right2', 'f4'])
            gendstats = pd.DataFrame(gendstats, columns = ['dscore', 'f4', 'Z'])
            self.detail = {'details': details, 'dscores': dscores, 'gendstats': gendstats}
            
    def __init__(self, logfile):
        
        qplog.__init__(self, logfile)
        self._get_anchor_line()
        if 'end_line' not in self._anchor_line.keys():
            raise ValueError('No EOF found in {}. Check if qpAdm was properly run.'.format(self.logfile))
        
        #self._parse_param(self._anchor_line['param_start_line'], self._anchor_line['version_line'])
        #self._parse_version(self._anchor_line['version_line'])
        #if self.version is None:
        #    raise ValueError('Please check if this is a qpAdm log.')
        
        self.left_pops = self._parse_left_right_pops(self._anchor_line['leftpop_start_line'] + 1)
        #self.right_pops = self._parse_left_right_pops(self._anchor_line['rightpop_start_line'] + 1)
        #self._parse_sample_size(self._anchor_line['sample_size_start_line'], self._anchor_line['data_info_start_line'])
        self._parse_coverage(self._anchor_line['coverage_start_line'] + 1)
        #self._parse_data_info(self._anchor_line['data_info_start_line'], self._anchor_line['coverage_start_line'])
        self._parse_data_info(self._anchor_line['used_data_info_start_line'], self._anchor_line['used_data_info_end_line'] + 1)
        
        if len(self.left_pops) == 2:
            self._parse_summary(self._anchor_line['summary_start_line'], self._anchor_line['summary_start_line'] + 11)
        else:
            self._parse_summary(self._anchor_line['summary_start_line'], self._anchor_line['summary_end_line'])
        
        #self._parse_details(self._anchor_line['detail_start_line'], self._anchor_line['end_line'])

### Parsing all tests

# Define relevant populations
a15 = ["Mota", "Ust_Ishim", "Kostenki14", "GoyetQ116-1", 
       "Vestonice16", "MA1", "ElMiron", "Villabruna", 
       "EHG", "CHG", "Iran-N", "Natufian", "Levant_N", 
       "WHG", "Anatolia_N"]

a15plus = a15 + ["Sar-Nur"]

post_Nur_inds = ["MSR002", "MSR003",
                 "VIL004", 'VIL006', 'VIL007', 'VIL009', 'VIL010', "VIL011",
                 "AMC001", "AMC005", "AMC014",
                 "COR001", "COR002",
                 "SNN001", "SNN002", "SNN003", "SNN004"]

mod_Sar = ["Cag", "Car", "Cam", "Ori", "Ogl", "Nuo", "Sas", "Olb"]
mod_Sar_inds = [p + str(i) for p in mod_Sar for i in range(1,6)]

sources_dict = {'Sard.Nuragic': ['Sar-Nur'],
           'N.African': ['Morocco_LN', 'Morocco_EN', 'Guanche'],
           'E.Mediterranean': ['Levant-BA', 'Iran-CA', 'Canaanite',
                               'Anatolia-BA', 'Minoan-BA', 'Myc-BA'],
           'Beaker': ['Italy_North-Bk', 'CE-EBA', 'CE-LBA',
                      'Iberia-BA', 'Balkans-BA'],
           'Sicily.BronzeAge': ['Sicily-Bk']}

pools = sources_dict.keys()


# utility functions for parsing
def parse_log_file(path):

    try :
        qpAdmlog = qpAdm_log(path)
        
    except ValueError:
        return [path, 'No EOF line']
    
    except FileNotFoundError:
        return [path, 'No file']
    
    else:
        left_pops = qpAdmlog.left_pops
        n = len(left_pops)
        sd = qpAdmlog.summary['std. errors']
        coverage = qpAdmlog.coverage['coverage'].values[:n]
        if n > 2:
            full_model = qpAdmlog.summary['model_comparison'].loc[[0]]
            full_model = full_model.values.tolist()[0][2:]
            return [path] + left_pops + list(coverage) + full_model[:-1]  + sd + [full_model[-1]]
        else:
            summ = qpAdmlog.summary['summ']
            return [path] + left_pops + list(coverage) + summ[-3:]

def produce_pool_combination(n, pools = pools):
    return  list(it.combinations(pools, n))

def parse_n_way_model(n, target_list = post_Nur_inds + mod_Sar_inds, sources_dict = sources_dict):
    '''
    Parse all n way models
    '''
    pool_combination = produce_pool_combination(n)
    
    results = []
    count = 0
    
    print('parsing {}-way model ... '.format(n))

    for c in pool_combination:
        if 'Sard.Nuragic' in c:
            out_dir = "NuragicSource/"
            rightpops = a15
            sources_list = [sources_dict[c[i]] for i in range(n)]
            sources_list = list(it.product(*sources_list))
            for t in target_list:
                for s in sources_list:
                    leftpops = [t] + list(s)
                    output_file = out_dir + '.'.join(leftpops) + ".log"
                    results.append(parse_log_file(output_file))
                    if results[-1][-1] == 'No EOF line':
                        print(results[-1][0])
                    count += 1
                    
        else:
            out_dir = "NuragicOutgroup/"
            rightpops = a15plus
            sources_list = [sources_dict[c[i]] for i in range(n)]
            sources_list = list(it.product(*sources_list))
            for t in target_list:    
                for s in sources_list:
                    leftpops = [t] + list(s)
                    output_file = out_dir + '.'.join(leftpops) + ".log"
                    results.append(parse_log_file(output_file))
                    count += 1
                    
    results_partial = [r for r in results if r[1] == 'No EOF line']
    results_miss = [r for r in results if r[1] == 'No file']
    results = [r for r in results if r[1] not in ['No EOF line', 'No file']]
    
    print('parsed {} files.'.format(count))
    print('{1}/{0} succeeded; {2}/{0} not finished; {3}/{0} not run.\n'.format(count, len(results), len(results_partial), len(results_miss)))
    
    return results, results_partial, results_miss


# save results into a error log file and a table
#if os.path.isfile('proximal_not_finished_jobs.txt'):
#    os.remove('proximal_not_finished_jobs.txt')

    

for i in range(1,2):
    succ, part, miss = parse_n_way_model(i)
    with open('proximal_ancinds/proximal_not_finished_jobs.txt', 'a') as file:
        for p in part:
            file.write('\t'.join(p) + '\n')
        for m in miss:
            file.write('\t'.join(m) + '\n')
         
    with open('proximal_ancinds/proximal_{}way_model.csv'.format(i), 'w') as file:
        if i == 1:
            header = ['file', 'target', 'source1', 'cov_target', 'cov1', 'pval', 'coef1', 'sd1']
        else:
            header = ['file', 'target'] + ['source'+ str(i) for i in range(1,i+1)] + ['cov_target']
            header += ['cov'+ str(i) for i in range(1,i+1)]
            header += ['df', 'chisq', 'pval']
            header += ['coef'+ str(i) for i in range(1,i+1)]
            header += ['sd'+ str(i) for i in range(1,i+1)]
            header += ['status']
        
        file.write(','.join(header) + '\n')
        for su in succ:
            file.write(','.join([str(s) for s in su]) + '\n')

