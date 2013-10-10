#!/usr/bin/env python 

import sys
from pprint import pprint

def parse_header(filename):
    modelinfo = {'pars':[]}
    pardesc = None
    gather_text = None
    with open(filename,"rt") as fid:
        for line in fid:

            # If line contains "<text>" then join lines together for
            # until line also contains "</text>"
            if gather_text:
                gather_text.append(line)
                if "</text>" in line:
                    line = "".join(gather_text)
                    gather_text = None
                else:
                    continue
            if "<text>" in line and not "</text>" in line:
                gather_text = [line]
                continue

            if "[DISP_PARAMS]" in line:
                modelinfo["disperse"] = getpars(line, ",")
            elif "[PYTHONCLASS]" in line:
                modelinfo["name"] = line.split("=")[1].strip().replace('Model','')
            elif "///" in line:
                pardesc = line.split("///")[1].strip()
                if pardesc.endswith("]"):
                    pardesc = pardesc.split("[",1)[0].strip()
            elif "[DEFAULT]" in line:
                _,parname,value_unit = line.split('=')
                parname = parname.strip()
                value_unit = value_unit.split()
                parvalue = float(value_unit[0])
                parunit = value_unit[1] if len(value_unit)>1 else ""
                if parunit.startswith("[") and parunit.endswith("]"):
                    parunit = parunit[1:-1]

                if not pardesc: pardesc = parname
                modelinfo["pars"].append({
                    'name': parname, 
                    'default': parvalue, 
                    'units': parunit, 
                    'description': pardesc,
                    'min': "0",
                    'max': "+DBL_INF"})
                pardesc = None
            elif "[FIXED]" in line:
                modelinfo["fixed"] = getpars(line)
            elif "[ORIENTATION_PARAMS]" in line:
                modelinfo["orient"] = getpars(line)
            elif "[MAGNETIC_PARAMS]" in line:
                modelinfo["magnetic"] = getpars(line)
            elif "[DESCRIPTION]" in line:
                line = line.replace('//','')
                modelinfo["description"] = trimdocs(gettext(line))
                #print "description:",cmultiline(modelinfo["description"])

    # Check that all parameters exist
    pars = set(p['name'] for p in modelinfo['pars'])
    for k in 'fixed','orient','magnetic','disperse':
        if set(modelinfo.get(k,[])) - pars:
            print pars, modelinfo[k]
            raise ValueError("missing parameters from %s: "%filename)
    return modelinfo

def gettext(line):
    line = line.split("=",1)[1].replace("<text>","").replace("</text>","")
    return line

def getpars(line, sep=";"):
    line = line.replace('*','')
    words = [w.strip() for w in gettext(line).split(sep)]
    return [w for w in words if w and not w.endswith('.width')]

def trimdocs(docstring):
    if not docstring:
        return ''
    lines = docstring.expandtabs().splitlines()

    # Determine minimum indentation (first line doesn't count):
    indent = sys.maxint
    for line in lines[1:]:
        stripped = line.lstrip()
        if stripped:
            indent = min(indent, len(line) - len(stripped))

    # Remove indentation (first line is special):
    trimmed = [lines[0].strip()]
    if indent < sys.maxint:
        for line in lines[1:]:
            trimmed.append(line[indent:].rstrip())

    # Strip off trailing and leading blank lines:
    while trimmed and not trimmed[-1]:
        trimmed.pop()
    while trimmed and not trimmed[0]:
        trimmed.pop(0)
    return '\n'.join(trimmed)

def cmultiline(string):
    return '\\\n' + "\\n\\\n".join(string.splitlines()) + '\\\n';

def gencdef(modelinfo):
    print '#include <ModelInfo.h>'
    print ''
    print 'ParameterInfo param_infos[] = {'
    print '   // name, description, units, default, min, max, flags'
    for p in modelinfo['pars']:
        flags = []
        if p['name'] in modelinfo.get('magnetic',[]):
            flags.append('PF_Magnetic')
        if p['name'] in modelinfo.get('orient',[]):
            flags.append('PF_Orientation')
        if p['name'] in modelinfo.get('disperse',[]):
            flags.append('PF_Polydisperse')
        if not flags: flags.append('PF_None')
        p = p.copy()
        p['flags'] = '|'.join(sorted(flags))
        print '   { "%(name)s", "%(description)s", "%(units)s", %(default)s, %(min)s, %(max)s, %(flags)s },'%p
    print '};'
    print ''
    print 'ModelInfo model_info('
    print '    "%s",'%modelinfo['name']
    print '    "%s",'%cmultiline(modelinfo['description'])
    print '    GetParameterCount(param_infos),'
    print '    param_infos);'

def gendisperser(modelinfo):
    print ''
    print '#include <disperser.h>'
    print 'typedef struct {'
    for p in modelinfo['pars']:
        print '  double %s;'%p['name']
    print '} Parameters;'
    print """
class Model: public Disperser {
public:
	Model(int pindex[], Weight pars[]) : _endpts(pindex), _pars(pars) {}

	double 
	formV(double dp[]) {
		Parameters *p = (Parameters *)dp;
		return
	}

	double
	formQ(double dp[], double q) {
		Parameters *p = (Parameters *)dp;
		return 
	}

	double
	formQxy(double dp[], double qx, double qy) {
		Parameters *p = (Parameters *)dp;
		return formQ(sqrt(qx*qx+qy*qy));
	}

	double
	formQxyz(double dp[], double qx, double qy, double qz) {
		Parameters *p = (Parameters *)dp;
		return
	}

	double
	formER(double dp[]) {
		Parameters *p = (Parameters *)dp;
		return
	}

	double
	formVR(double dp[]) {
		Parameters *p = (Parameters *)dp;
		return 1.0;
	}

} ;

//CExport void* get_model_info() { return NULL; }
//CExport void* create_model(void* data) { return NULL; }
//CExport void destroy_model(void* ptr) {}
CExport void calculate_q(int pindex[], Weight pars[], size_t nq, double iq[], double q[]) {
	Model(pindex, pars).calc_q(nq, q, iq);
}
CExport void calculate_qxqy(int pindex[], Weight pars[], size_t nq, double iq[], double qx[], double qy[]) {
	Model(pindex, pars).calc_qxqy(nq, qx, qy, iq);
}
CExport void calculate_qxqyqz(int pindex[], Weight pars[], size_t nq, double iq[], double qx[], double qy[], double qz[]) {
	Model(pindex, pars).calc_qxqyqz(nq, qx, qy, qz, iq);
}
CExport double calculate_ER(int pindex[], Weight pars[]) {
	Model(pindex, pars).calc_ER();
}
CExport double calculate_VR(int pindex[], Weight pars[]) {
	Model(pindex, pars).calc_VR();
}
"""
    

def main():
    for f in sys.argv[1:]:
        print "===== file:",f
        #pprint(parse_header(f)) 
        info = parse_header(f)
        gencdef(info)
        gendisperser(info)

if __name__ == "__main__":
    main()
