"""
c2pap-web-portal batch jobs, by Antonio Ragagnin (ragagnin@lrz.de), 2018.

many thanks to Leonard Bauer for his beta testing.

Mini Documentation
==================

Once you have a list of clusters is your dataset.csv, you can run a batch of PHOX jobs via
 
python c2pap_batch.py -u 'your@email' -f dataset.csv -s PHOX -p mode='ICM only'  instrument='eROSITA (A=1000 FoV=60)' t_obs_input=1000 img_z_size=200 simulate=1

or to run SMAC:

python c2pap_batch.py -f dataset.csv -u 'your@rmail'  -s SMAC -p content="bolometric x-ray luminosity" IMG_SIZE=512 IMG_Z_SIZE=5000 PROJECT='along y, xz plane' 

To save computing time (and your time) it is possible to tell the script to recycle your pre-existing jobs with the same parameters (very useful when you run batch scripts that intersects results from old scripts). For this purpose, add the parameters -e -j <name_of_a_cache_file> e.g.
 
python c2pap_batch.py -u 'your@email' -f dataset.csv -s PHOX -p mode='ICM only'  instrument='eROSITA (A=1000 FoV=60)' t_obs_input=1000 img_z_size=200 simulate=1 -e -j jobs_cache.pk
 
The first time you execute it, it will take a lot of time because it will have to load all of your jobs. The next times it will be very fast (it will load them from the file job_cache.pk).

"""
import getpass
import mechanize
import argparse
import sys
import csv
import time
import os
import re 
import pickle 
import signal
import time
import json
import urllib


BASE = "https://c2papcosmosim.uc.lrz.de" 
__version__= "1.0 [22 november 2018]"
__author__= "Antonio Ragagnin"
globy = {}

        
def find_between( s, first, last ):
    "just return the first sub-string of `s` between string delimeters `first` and `last`"
    start = s.index( first ) + len( first )
    end = s.index( last, start )
    return s[start:end]


def find_between_r( s, first, last ):
    "as find_between, but return the last sub-string"
    start = s.rindex( first ) + len( first )
    end = s.rindex( last, start )
    return s[start:end]



def get_between_tags(res, tag, flags=None):
    "given a HTML content in `res`, returns a list of string for each sub-string of `res` between the HTML tags <`tag`> and</`tag`> "
    len_tag = len(tag)
    if flags is None:
        return [res[x.start()+len_tag+2:x.end()-(len_tag+3)] for x in re.finditer('<%s>.*?</%s>'%(tag,tag),res)]
    else:
        return [res[x.start()+len_tag+2:x.end()-(len_tag+3)] for x in re.finditer('<%s>.*?</%s>'%(tag,tag),res, flags=flags)]


def get_jobs_ids(br):
    "open the url with the job lists and parse the HTML to find all job ids"
    res = br.open(BASE+"/jobs").read()
    #the job ids are inside a <code></code> HTML tag
    jobids = get_between_tags(res, 'code')
    no_error_jobids = []
    for jobid in jobids:
        if find_between(res, jobid,'</tr>' ).find("ERROR")==-1:
            no_error_jobids.append(jobid)
    return no_error_jobids


def login(br, username, password):
    "login to the web portal"
    br.open(BASE+"/login")
    br.select_form(nr=0)
    br.form.find_control("email").value = username
    br.form.find_control("password").value = password
    res = br.submit().read()
    if 'password' in res:
        raise Exception("Login failed for user=%s."%(username))


def get_job_param(res):
    """`res` should contain the HTML of the web page of a job. 
    The following function parses its table with the job parameter names and value
    (e.g. instrument = eRosita) and returns it as a dictionary d[job_parameter]=job_value"""
    d={}
    trs= get_between_tags(res, 'tr', flags=re.DOTALL)
    for tr in trs:
        tds= get_between_tags(tr, 'td', flags=re.DOTALL)
        if len(tds)==2:
            d[tds[0]]=tds[1]
    return d
  

def printf(format,err=False):
    fd = sys.stderr if err else sys.stdout
    fd.write(format)

def log(*s,**l):
    for x in s:
        sys.stderr.write(str(x)+" ")
    for x in l.keys():
        sys.stderr.write(str(x)+"="+str(l[x]))
    sys.stderr.write("\n")



def download(br, naming, jobid, job):
    """ Given the `jobid` and the dictionary of the job parameters `job` 
    (`job` is a dictionary computed by `get_job_param`),
    it connects to the FTP and download job results. 
    The folder name is decided using the template string `naming`.
    For instance if `naming` is taken by the command line parameter -n,
    it will have a default value of
    '{simulation_name}/{SnapNum}/{cluster_id}/{application_name}/{jobid}'  """
    data = dict(job)
    data["jobid"]=jobid
    name = naming.format(**data)
    os.system("mkdir -p %s"%(name))
    code_to_download = "wget ftp://129.187.239.195:/%s/%s.tar.gz -qO -|  tar -xz -C %s/ --strip-components=1 %s"%(jobid, jobid, name, jobid)
    log("Executing: %s"%code_to_download)
    os.system(code_to_download)
    with open("%s/job.par"%(name),"w") as f:
        for k in data.keys():
            f.write("%s %s\n"%(k, data[k]))
    log("Written: %s"%name)

def wait(br, jobid):
    "given a `jobid`, this function check the list of jobs until `jobid` results as 'COMPLETED' "
    res = br.open(BASE+"/jobs").read()
    status=None
    while status!="COMPLETED":
        row = find_between(res, '<code>%s</code>'%(jobid), '</tr>')
        status = find_between(row, '<td>', '</td>')
        log("status: ",status)
        if status!="COMPLETED":
            log("now waiting %ds..."%(globy["time_sleep"]))
            time.sleep(globy["time_sleep"])
            res = br.open(BASE+"/jobs").read()

    return jobid

def fill(br,f,v, debug=False, byval=False):
    """very tricky and foggy function. 
    Given a mechanize browser `br`, it tries to set the value `v` 
    in the field <form> field named `f`.
    There are several switches in this function because addigment 
    is different wehn you set a input-text, a check-box, a combo-box and so on """
    if br.form.find_control(f).type == 'select' and byval==False:
        item_set = False
        for item in br.form.find_control(f).items:
            if v == item.attrs['label']:
                br.form.find_control(f).value=[item.attrs['value']]
                print(f,"-> [",v,']=',item.attrs['value'])
                item_set = True
        if item_set == False:
            raise Exception("Error: field %s has no value %s"%(f,v))
    else:
        print(f,"->",v)
        br.form.find_control(f).readonly = False
        if br.form.find_control(f).type == 'select':
            br.form.find_control(f).value = [v]
        elif    br.form.find_control(f).type == 'checkbox':
            br.form.find_control(f).value = ['on'] if v[0] in ['t','T','y','Y','1','o'] else []
            #print(br.form.find_control(f).checked)
        else:
            br.form.find_control(f).value = v

def submit(br, cluster_id, args):
    """given a `cluster_id`, it submit a job with the values of the form set in `args`
    (e.g. args['content']='content="bolometric x-ray luminosity')."""
    log("Opening: %s"%BASE+"/map/%s"%(args.service.lower()))
    br.open(BASE+"/map/%s"%(args.service.lower()))
    br.select_form(nr=0)
    fill(br, "cluster_id", str(int(cluster_id)), debug=True)
    for k in args.ps.keys():
        if k=='redshift': continue
        v=args.ps[k]
        if v.is_option_number:
            fill(br, k, str(v), debug=True, byval=True)
        else:
            fill(br, k, str(v), debug=True, byval=False)
    #sys.exit(0)
    res = br.submit(name="submit_field").read()
    log("Job submitted.")
    return  find_between(res, "ID <a href='/jobs/",'/')

def cmpf(a,b,precision=0, debug=False):
    "compare two floats within a precision. I use it to check if two jobs have the same parameters."
    s='%%.%df'%(precision)
    if debug:
        pass
        print( s%float(a) == s%float(b), s%float(a) , s%float(b) )
    return  s%float(a) == s%float(b)

def compare(job, args, debug=True, service_data=None):
    """check if the job parameters that were set as input to this program match
    the job-parameter of an already exiting job `job`. As you can imagine 
    this function is very delicate as for instance the user send a job with
    comoving distances, but its job-parameter  is further stored in physical units.
    So there are several conversions.   """
    if args.service=="SimCut":
        
        r500factor1 = args.ps['r500factor']
        z_image_size1 = args.ps['IMG_Z_SIZE']
        r500factor2 = float(job['IMG_XY_SIZE'])/ float(job['R_VIR'])/2.
        z_image_size2 = job['IMG_Z_SIZE']
        #print(z_image_size1,z_image_size2,r500factor1, r500factor2)
        return  cmpf(z_image_size1,z_image_size2, debug=debug) and cmpf(r500factor1, r500factor2, precision=1, debug=debug)
    if args.service=="SMAC":
        c = True
        c = c and (args.ps['content']==job['content'])
        c = c and cmpf(args.ps['IMG_Z_SIZE'], job['IMG_Z_SIZE'], debug=debug)
        project = None

        if args.ps['PROJECT'].is_option_number:
            project = str(service_data['PROJECT'][int(args.ps['PROJECT'])]['value'])
        else:
            for p in service_data['PROJECT']:
                if p['name']==args.ps['PROJECT']:
                    project = str(p['value'])
                    break
        if project is None:
            raise Exception("Unable to find the value of PROJECT='%s'"%(args.ps['PROJECT']))
            
        c = c and project == job['PROJECT']
        c = c and cmpf( float(args.ps['r500factor'])/(1.+float(args.ps['redshift']))/(float(job['HUBBLE'])/100), float(job['IMG_XY_SIZE'])/ float(job['R_VIR'])/2., debug=debug)
        return c

    if args.service=="PHOX":
        c = True
        mode_id = ['ICM only','AGN only','ICM+AGN'].index(args.ps['mode'])+1
        c = c and (str(mode_id)==job['mode'])
        c = c and cmpf(args.ps['img_z_size'], job['img_z_size'], debug=debug)
        c = c and cmpf(args.ps['t_obs_input'], job['T_obs'], debug=debug)
        c = c and str(args.ps['instrument'].split(' ')[0]==job['instrument'])
        c = c and str(args.ps['simulate'][0]  in ['t','T','y','Y','1','o'] and job['simulate']=='1')

        return c
    else:
        raise ValueError



def save_cache(cache_jobs, cache):
    "save all jobs data in a cache file, so next time we load it faster"
    #we disable the CTRL+C signal while saving, we do not want a corrupted file
    s = signal.signal(signal.SIGINT, signal.SIG_IGN)
    with open( cache_jobs, "wb" )  as f:
        pickle.dump( cache, f )
    #re-enable CTRL+C
    signal.signal(signal.SIGINT, s)


def query(br, query, page, limit):
    "personal use"
    j = br.open(BASE+"/query/raw?"+urllib.urlencode((('q',query),('p',page),('l',limit)))).read()
    o = json.loads(j)
    return o

def get_cache(cache_jobs):
    "load file from cache if exists, otherwise create one." 
    try:
        print("open %s"%(cache_jobs))
        with open( cache_jobs, "rb" )  as f:
            pass
    except:
        cache =  {}
        print("write %s"%(cache_jobs))
        save_cache(cache_jobs, cache)
        return cache

    print("read %s"%(cache_jobs))
    with open( cache_jobs, "rb" )  as f:
        cache = pickle.load(f)
    return cache


gcache={"x":{}}
def  get_job_or_cache(br, jobid, cache_jobs):
    """given a job with `jobid`, either load it from a cahce file or,
    if not there, load it from the web. Then the job is saved in the cache"""
    if cache_jobs is None:
        log("Loading: %s"%(jobid))
        res = br.open(BASE+"/jobs/%s/show"%(jobid)).read()
        p = get_job_param(res)
        return p

    if "jobs" not in  gcache["x"]:
        gcache["x"] = get_cache(cache_jobs)
    if "jobs" not in  gcache["x"]:
        gcache["x"]["jobs"] = {}

    cache = gcache["x"]
    jobs = cache["jobs"]
    if jobid not in jobs:
        log("Loading: %s"%(jobid))
        res = br.open(BASE+"/jobs/%s/show"%(jobid)).read()
        p = get_job_param(res)
        jobs[jobid]=p
        save_cache(cache_jobs, cache)
    return jobs[jobid]

class Parameter(str):
    "store job parameters"
    def __new__(cls, value, *args, **kwargs):
        return super(Parameter, cls).__new__(cls, value)
    def __init__(self, value, is_option_number=False):
        self.is_option_number=is_option_number

def get_sims_and_snaps(br, snap_id):
    """given the snap_id inside dataset.csv, we parse the web-portal HTML
    in order to find which simulation we are using.
    Then, the browser select the simulation and the snap, so the jobs
    are executed there."""
    res = br.open(BASE+'/map/find').read()
    j = json.loads(find_between(res, "myApp.constant('angulardata',", ")\n"))
    simulations = j["simulations"]
    print (simulations)

    select_simulation = None
    select_snap = None
    redshift = None
    #print json.dumps(j, indent=4, sort_keys=True)

    for simulation in simulations:
        br.select_form(nr=0)
        #log("test %s with id=%d"%(simulation["name"], simulation["id"]))
        fill(br, "simulation", str(simulation["id"]), byval=True)
        res = br.submit().read()
        j = json.loads(find_between(res, "myApp.constant('angulardata',", ")\n"))
        snaps = j["snaps"]
        for snap in snaps:
            if snap["id"] == snap_id:
                br.select_form(nr=0)
                fill(br, "snapshot", str(int(snap_id)), byval=True)
                select_snap = snap["name"]
                select_simulation = simulation["name"]
                redshift = float(snap["redshift"])
                log("Found %s/%s"%(select_simulation,select_snap)) 
                br.submit()
                #return str(simulation["id"]), str(int(snap_id))
                return  str(simulation["id"]), str(int(snap_id)), select_simulation, select_snap, redshift
    raise Exception("Unable to set simulation/snap")

class MyBrowser(mechanize.Browser):
    """Extend mechanize browser to keep retrying when we get a 502 error.
    I do this because the web portal sometimes give a 502 error and if
    if happens while you are processing thousand of jobs, it is annoying
    to wake up in the morning and discover that the code crashed because of a 502 error."""
    def __new__(cls, value, *args, **kwargs):
        return super(MyBrowser, cls).__new__(cls, value)
    def open(br, page):
        printf("Opening %s\n"%(page),err=True)
        while True:
            try:
                
                return mechanize.Browser.open(br, page)
            except Exception as e:
                if e.code == 502:
                    printf("Bad Gateway :( retrying in 60 seconds..\n",err=True)
                    time.sleep(globy["time_sleep"])
                else:
                    raise
def main():
    log ("")
    log ("===================================================")
    log ("=                                                 =")
    log ("= c2pap-web-portal batch jobs                     =")
    log ("= By Antonio Ragagnin, 2018.                      =")
    log ("= questions: ragagnin@lrz.de/ragagnin@usm.lmu.de  =")
    log ("=                                                 =")
    log ("===================================================")
    log ("")
    log ("Use option -h for help.")
    log ("")
    log ("version: %s."%(__version__))
    log ("")
 
    from argparse import RawTextHelpFormatter
    parser = argparse.ArgumentParser(description='c2pap webportal job batches', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-f', '--file', type=str,  help='"dataset.csv" returned by webportal', default='dataset.csv')
    parser.add_argument('-u', '--username', type=str,  help='If not used, the email will be read from standard input.', required=True)
    parser.add_argument('-x', '--password', type=str,  help='If not used, the password will be safely read from standard input.')
    parser.add_argument('-v', '--verbose',  action="count", help='Display HTTPs requests')
    parser.add_argument('-w', '--weak-ssl', action='store_true', help="Doesn't check SSL certificate.", default=True)
    parser.add_argument('-e', '--existing-jobs', action='store_true',  help="Search for jobs with same parameters. If found, download its data.", default=False)
    parser.add_argument('-s', '--service', choices=['SMAC','SimCut','PHOX','query'], help="Choose a service between SMAC, SimCut, PHOX. Check the lower/upper case.", required=True)
    parser.add_argument('-t', '--time-sleep', help="How frequently check the status of the job in seconds. Default=60", default=60, type=int)
    parser.add_argument('-a', '--auto-download', help="Download data after execution, default=True", default=True, type=bool)
    parser.add_argument('-p', '--params', help="""Form parameters.
    SimCut:
    IMG_Z_SIZE
    r500factor

    SMAC:
    content
    IMG_SIZE
    IMG_Z_SIZE
    PROJECT
    r500factor

    PHOX:
    mode ['ICM only','AGN only','ICM+AGN'] 
    instrument (-1 for generic) or 'eROSITA (A=1000 FoV=60)' ..
    instrument_a (only if generic)
    instrument_fov (only if generic)
    t_obs_input e.g. 1000
    img_z_size e.g. 2000
    simulate use '1' or '0'

    query:
    query
    page
    limit
    """,  nargs='+', default=["IMG_Z_SIZE=200","r500factor=1.0"])
    parser.add_argument('-j', '--cache-jobs', help="Cache jobs values in file", default=None, type=str)# "jobs.cache", type=str)
    parser.add_argument('-n', '--naming', help="""
    Where to save data. Use interpolation with %%(variable name). 
    Variable name can be  %%(jobid) or any parameter on https://c2papcosmosim.uc.lrz.de/jobs/%%(jobid)/show.

    default: {simulation_name}/{SnapNum}/{cluster_id}/{application_name}/{jobid}
    
""", default="{simulation_name}/{SnapNum}/{cluster_id}/{application_name}/{jobid}", type=str)

    args = parser.parse_args()

    globy["time_sleep"] = args.time_sleep
    globy["auto_download"] = args.auto_download
    args.ps = {}
    for p in args.params:
        try:
            if ':=' in p:
                raise Exception("assignment with ':=' is no more possible!")
                k,v = p.split(':=', 1)
                args.ps[k]=Parameter(v, is_option_number=True)
            elif '=' in p:
                k,v = p.split('=', 1)
                args.ps[k]=Parameter(v)


        except:

            printf("Error evaluating %s\n"%(p),err=True)
            sys.exit(1)

    if args.weak_ssl:
        import ssl
        ssl._create_default_https_context = ssl._create_unverified_context

    log ("Initializing browser...")
    br =MyBrowser()

    br.set_handle_equiv(True)
    br.set_handle_gzip(True)
    br.set_handle_redirect(True)
    br.set_handle_referer(True)
    br.set_handle_robots(False)
    if args.verbose:
        if args.verbose>=1:
            br.set_debug_http(True)
        if args.verbose>=2:
            br.set_debug_redirects(True)
        if args.verbose>=3:
            br.set_debug_responses(True)

    br.set_handle_robots(False)
    br.addheaders = [('User-agent', 'Mozilla/5.0')] #we are liars mwhahahah

                  

    log ("Logging in...")
    if args.username == None:
        args.username = raw_input('username: ')
    if args.password == None:
        args.password = getpass.getpass('password: ') #safely ask for password if not enterd via command line.
    login(br, args.username, args.password)
    args.password = "" #after login, we  remove password from memory
    del args.password

    
    service_data = {}

    if args.service=="SMAC":
        import yaml
        service_data = yaml.load(br.open(BASE+"/static/smac.yml"))

    if args.service=="PHOX":
        import yaml
        service_data = yaml.load(br.open(BASE+"/static/phox.yml"))


    if args.service == "query":
        if 'query' not in args.ps:
            log("specify --param query='query'")
            sys.exit(1)
        limit = args.ps['limit'] if 'limit' in args.ps else 100 
        page =  args.ps['page'] if 'page' in args.ps else 0 
        j = query(br, args.ps['query'], page, limit)
        if 'error' in j:
            printf(j['error'],err=True)
            printf("\n",err=True)
            sys.exit(1)
                
        if 'header' not in j:
            printf(json.dumps(j),err=True)
            printf("\n",err=True)

            sys.exit(1)
        h = j['header']
        rs = j['rows']
        print '#',
        for label in h:
            print label,
        print
        for row in rs:
            for cell in row:
                print cell,
            print
        sys.exit(0)
    jobs = {}
    if args.cache_jobs:
        br.cache = None


    log ("Loading existing jobs...")
    if args.existing_jobs:
        job_ids = get_jobs_ids(br)
        for jobid in job_ids:
            job = get_job_or_cache(br, jobid, args.cache_jobs)
            jobs[jobid] = job
    

    log ("Reading clusters list from: '%s' ..."%(args.file))
    clusters=[]
    with open(args.file, 'rb') as f:
        clusters = [{k: float(v) for k, v in row.items()}
             for row in csv.DictReader(f, skipinitialspace=True)] #no idea, I just found this code on stack overflow 
    log ("Read %d clusters."%(len(clusters)))


    simulationid, snapid, simulation_name, snap_name, redshift =  get_sims_and_snaps(br,clusters[0]["snap_id"])
    #print(simulation_name, snap_name,  redshift)
    args.ps["simulation"]=Parameter(simulationid, is_option_number=True)
    args.ps["snapshot"]=Parameter(snapid, is_option_number=True)
    args.ps["redshift"]=Parameter(redshift, is_option_number=False)


    for cluster in clusters:
            found_job = None
            cluster_id = cluster["uid"]
            log("")
            log("Cluster: %s/%s/%d, z=%.1f internalid:%d"%(simulation_name, snap_name, cluster_id,redshift, cluster['# id']))

            if args.existing_jobs:
                for jobid in jobs.keys():
                    job = jobs[jobid]
                    if "cluster_id" not in job:
                        continue
                    #print(int(cluster['# id']),simulation_name,job["simulation_name"], snap_name, job["SnapNum"] , cluster_id,  job["cluster_id"] ,job["application_name"], args.service.lower(),
                    #snap_name==job["SnapNum"] ,job["simulation_name"] == simulation_name, cmpf(job["cluster_id"], str(cluster_id)),job["application_name"]==args.service.lower())
                    #print(jobid, snap_name, job["SnapNum"],  job["simulation_name"] , simulation_name, job["application_name"], args.service.lower())
                    if (snap_name==job["SnapNum"] and 
                    job["simulation_name"] == simulation_name and 
                        cmpf(job["cluster_id"] , str(cluster_id) )and  
                    job["application_name"]==args.service.lower()):
                        #print("compate", job, args)
                        #print("FOUN?")
                        found = compare(job, args, service_data=service_data)#,debug=True)
                        if (found):
                            found_job = jobid
                        #print("found",found)
            if  found_job is None:
                #print cluster
                log("Preparing job on (internal) clusterid=%d center=[%.0f,%.0f,%.0f].."%(int(cluster['uid']), cluster["x"], cluster["y"], cluster["z"]))
                found_job = submit(br, cluster["# id"], args)
                log("New job: %s."%( found_job) )
            else:
                log("Old job: %s."% (found_job) )

            wait(br, found_job)
            download(br, args.naming, found_job,  get_job_or_cache(br, found_job, args.cache_jobs))
            #sys.exit(0)


    
    
    
if __name__=="__main__":
    main()
