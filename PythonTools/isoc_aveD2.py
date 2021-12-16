import numpy as np
import time


def main():
    time_start = time.time()
    dir_name = "isoc_dt0.1_data.T0.38.Equi.DPD"
    set_iffpath = "./{it_sample}/{dir_name}/D2/D2.{it_isoc}"
    set_offpath = "./{it_sample}/{dir_name}/D2/D2.ave"
    particle_number = 4320
    header = 9
    list_str_header = []
    isoc_num = 100
    sample_num = 20
    for it_sample in range(sample_num):
        list_ave_D2 = []
        for it_isoc in range(isoc_num):
            iffpath = set_iffpath.format(it_sample=it_sample + 1,
                                         dir_name=dir_name,
                                         it_isoc=it_isoc + 1)
            list_D2_isoc = D2_file_to_list(iffpath=iffpath, header=header)
            if len(list_D2_isoc) != particle_number:
                print(len(list_D2_isoc))
                raise Exception("wrong shape: it_list_D2")
            if it_isoc == 0:
                list_ave_D2 = list_D2_isoc
            else:
                for it_list_D2 in range(particle_number):
                    list_ave_D2[it_list_D2].strict_safe_add_D2(
                        list_D2_isoc[it_list_D2])
        for it_list_ave_D2 in range(particle_number):
            list_ave_D2[it_list_ave_D2].D2 /= isoc_num
            # print(list_ave_D2[it_list_ave_D2].D2)
        offpath = set_offpath.format(it_sample=it_sample + 1,
                                     dir_name=dir_name)
        with open(
                set_iffpath.format(it_sample=it_sample + 1,
                                   dir_name=dir_name,
                                   it_isoc=0 + 1), 'rt') as f:
            for i in range(header):
                list_str_header.append(f.readline())
            # print(len(list_str_header))

        with open(offpath, 'wt') as of:
            for i in range(header):
                print(list_str_header[i], file=of, end='')
            for i in range(particle_number):
                iline = list_ave_D2[i].to_string()
                # print(iline)
                print(iline, file=of)
        print('{} done!'.format(it_sample))
        time_end = time.time()
        print('time cost:', time_end - time_start, 's')


class Pa_D2:
    'particle with d2 information'
    pid = 0
    ptype = 0
    x = 0.
    y = 0.
    z = 0.
    ix = 0
    iy = 0
    iz = 0
    D2 = 0.

    def __init__(self, list=[0, 0, 0., 0., 0., 0, 0, 0, 0]):
        self.pid = list[0]
        self.ptype = list[1]
        self.x = list[2]
        self.y = list[3]
        self.z = list[4]
        self.ix = list[5]
        self.iy = list[6]
        self.iz = list[7]
        self.D2 = list[8]

    def add_D2(self, Pa_D2):
        self.D2 = self.D2 + Pa_D2.D2
        return self

    def strict_safe_add_D2(self, Pa_D2):
        safe_flag = (self.pid == Pa_D2.pid)\
                and (self.ptype == Pa_D2.ptype)\
                and (self.x == Pa_D2.x)\
                and (self.y == Pa_D2.y)\
                and (self.z == Pa_D2.z)\
                and (self.ix == Pa_D2.ix)\
                and (self.iy == Pa_D2.iy)\
                and (self.iz == Pa_D2.iz)
        if safe_flag == True:
            return self.add_D2(Pa_D2)
        else:
            raise Exception("particles not the same when add D2!")

    def idsafe_add_D2(self, Pa_D2):
        safe_flag = safe_flag = (self.pid == Pa_D2.id)
        if safe_flag:
            return self.add_D2(Pa_D2)

    def to_string(self):
        fm_ostr = "{:.0f} {:.0f} {:.6g} {:.6g} {:.6g} {:.0f} {:.0f} {:.0f} {:.6g}"
        ostr = fm_ostr.format(self.pid, self.ptype, self.x, self.y, self.z,
                              self.ix, self.iy, self.iz, self.D2)
        return ostr


def D2_file_to_list(iffpath, header):
    d2_it = np.loadtxt(iffpath, skiprows=header)
    total_row = d2_it.shape[0]
    total_col = d2_it.shape[1]
    list_D2 = []
    for it_row in range(total_row):
        itrow_pa_D2 = Pa_D2(d2_it[it_row][:])
        list_D2.append(itrow_pa_D2)
    return list_D2


if __name__ == "__main__":
    main()