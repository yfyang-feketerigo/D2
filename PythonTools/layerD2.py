from functools import total_ordering
import numpy as np
import time
import math


def main():
    time_start = time.time()
    dir_name = "isoc_dt0.1_data.T0.38.Equi.DPD"
    set_iffpath = "./{it_sample}/{dir_name}/D2/D2.ave"
    set_offpath = "./{it_sample}/{dir_name}/D2/D2.layer"

    layer_num = 20
    total_length = 19.31
    header = 9

    sample_num = 20
    for isample in range(sample_num):
        it_sample = isample + 1
        list_aveD2 = []
        iffpath = set_iffpath.format(it_sample=it_sample, dir_name=dir_name)
        list_aveD2 = D2_file_to_list(iffpath=iffpath, header=header)
        list_layerD2 = list_D2_to_list_layerD2(list_D2=list_aveD2,
                                               layer_num=layer_num,
                                               total_len=total_length)
        offpath = set_offpath.format(it_sample=it_sample, dir_name=dir_name)
        with open(offpath, 'wt') as of:
            print("layer D2 counter", file=of)
            for i in range(len(list_layerD2)):
                print(list_layerD2[i].D2info_to_string(), file=of)
        print('{} done'.format(it_sample))
    time_stop = time.time()
    print('done! time cost:', time_stop - time_start, 's')


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


class layerD2:
    layer_id = 0
    pa_counter = 0
    D2 = 0.
    total_layer_num = 0
    total_length = 0.

    def __init__(self, layer_id, D2, pa_coutner, total_layer_num,
                 total_length):
        self.layer_id = layer_id
        self.D2 = D2
        self.pa_counter = pa_coutner
        self.total_layer_num = total_layer_num
        self.total_length = total_length

    def D2info_to_string(self):
        fmt_str = "{:.0f} {:.6g} {:.0f}"
        return fmt_str.format(self.layer_id, self.D2, self.pa_counter)


def list_D2_to_list_layerD2(list_D2, layer_num, total_len):
    # COL_GRAD=3
    # COL_D2=8

    dlayer = total_len / layer_num
    list_layerD2 = []
    for i in range(layer_num):
        list_layerD2.append(
            layerD2(layer_id=0,
                    D2=0.,
                    pa_coutner=0,
                    total_layer_num=layer_num,
                    total_length=total_len))
    # print(list_layerD2[19].pa_counter)
    for i in range(len(list_D2)):
        layer = math.floor(list_D2[i].y / dlayer)
        list_layerD2[layer].D2 += list_D2[i].D2
        list_layerD2[layer].pa_counter += 1
        # print(list_layerD2[layer].pa_counter)
    for i in range(layer_num):
        list_layerD2[i].id = i
        list_layerD2[i].D2 /= list_layerD2[i].pa_counter
        # print(list_layerD2[i].pa_counter)
    return list_layerD2


if __name__ == "__main__":
    main()