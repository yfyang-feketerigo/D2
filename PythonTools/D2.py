import os

sample_start = 1
sample_stop = 20

isoc_start = 1
isoc_stop = 100

rcut = 2.5

isoc_dir = "isoc_dt1_wi50_data.T0.38.Equi.DPD"
pairstyle_t = "none"
pairstyle_t0 = "none"
boxtype_t = "tilt"
boxtype_t0 = "orthogonal"
ifname_t0 = "data.T0.38.Equi.DPD"
set_ifname_t = "{}.isco"
set_ifpath = "./{sampledir}/{isoc_dir}/"
set_ofname = "D2.{}"
for it_sample in range(sample_start, sample_stop):
    ifpath = set_ifpath.format(sampledir=it_sample, isoc_dir=isoc_dir)
    ofpath = ifpath + "D2/"
    for it_isoc in range(isoc_start, isoc_stop):
        ifname_t = set_ifname_t.format(it_isoc)
        ofname = set_ofname.format(it_isoc)

        set_run_D2 = "./D2.exe -r {r} -t {t} -0 {t0} -I {I} --boxtype_t {boxtype_t} --boxtype_t_ {boxtype_t_} -o {o} -O {O}"
        flag = os.system(
            set_run_D2.format(r=rcut,
                              t=ifname_t,
                              t0=ifname_t0,
                              I=ifpath,
                              boxtype_t=boxtype_t,
                              boxtype_t_=boxtype_t0,
                              o=ofname,
                              O=ofpath))
        if (0 != flag):
            print("erro code: {}".format(flag))
            exit(flag)
