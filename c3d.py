import ezc3d
import numpy as np
import os
import pandas as pd
from pathlib import Path
import shutil


file_address='input_your_file_address'

def get_file_paths_with_extension(folder_path, file_extension):
    # 列出文件夹中的所有文件和目录
    all_files = os.listdir(folder_path)

    # 过滤出特定后缀的文件
    specific_files = [f for f in all_files if f.endswith(file_extension)]

    # 生成文件的完整路径
    file_paths = [os.path.join(folder_path, f) for f in specific_files]

    return file_paths



folder_path = file_address
file_extension = '.c3d'
file_paths = get_file_paths_with_extension(folder_path, file_extension)
file_paths = [item.replace('\\', '/') for item in file_paths]
print(f'文件读取成功,共{len(file_paths)}个文件,正在处理中....')



for i in range(len(file_paths)):

    file=file_paths[i]
    c3d=ezc3d.c3d(file)

    markers=c3d['data']['points']
    num_frames=c3d['header']['points']['last_frame']-c3d['header']['points']['first_frame']+1
    num_markers=markers.shape[1]
    rate = c3d['parameters']['POINT']['RATE']['value'][0]
    labels = c3d['parameters']['POINT']['LABELS']['value']
    num_channel = len(labels)
    time = np.linspace(0, num_frames / rate, num_frames)
    rate = c3d['parameters']['ANALOG']['RATE']['value'][0]
    frame = c3d['header']['analogs']['last_frame'] - c3d['header']['analogs']['first_frame'] + 1




    with open(f'{i}.trc', 'w') as f:
        f.write("PathFileType\t4\t(X/Y/Z)\t" + file + "\n")
        f.write("DataRate\tCameraRate\tNumFrames\tNumMarkers\tUnits\tOrigDataRate\tOrigDataStartFrame\tOrigNumFrames\n")
        f.write(f"{rate:.2f}\t{rate:.2f}\t{num_frames}\t{num_markers}\tmm\t{rate:.2f}\t1\t{num_frames}\n")


        f.write("Frame#\tTime\t" + "\t".join([label + "\t\t" for label in labels]) + "\n")
        f.write("\t\t" + "\t".join(
            ["X" + str(i + 1) + "\tY" + str(i + 1) + "\tZ" + str(i + 1) for i in range(num_markers)]) + "\n")

        for frame_idx in range(num_frames):
            line_data = [f"{frame_idx + 1}", f"{time[frame_idx]:.5f}"]
            for marker_idx in range(num_markers):
                x = markers[0, marker_idx, frame_idx]  # X 轴保持不变
                y = markers[2, marker_idx, frame_idx]  # 新的 Y 取自 Z
                z = -markers[1, marker_idx, frame_idx]  # 新的 Z 取自负 Y
                line_data.extend([f"{x:.5f}", f"{y:.5f}", f"{z:.5f}"])  # 写入每个标记的X, Y, Z坐标
            f.write("\t".join(line_data) + "\n")
        f.close()


for i in range(len(file_paths)):

    file=file_paths[i]
    c3d=ezc3d.c3d(file)
    rate = c3d['parameters']['ANALOG']['RATE']['value'][0]
    frame = c3d['header']['analogs']['last_frame'] - c3d['header']['analogs']['first_frame'] + 1
    time = np.linspace(0, frame / rate, frame)
    channel = len(c3d['parameters']['ANALOG']['LABELS']['value'])
    force = c3d['data']['analogs']
    force_type=c3d['parameters']['FORCE_PLATFORM']['TYPE']['value'][0]


    if force_type == 2:
        num_force = int(channel /6)
        labels = []
        for j in range(num_force):
            a = f'{j}_ground_force_vx'
            b = f'{j}_ground_force_vy'
            c = f'{j}_ground_force_vz'
            d = f'{j}_ground_force_px'
            e = f'{j}_ground_force_py'
            f = f'{j}_ground_force_pz'
            g = f'{j}_ground_torque_x'
            h = f'{j}_ground_torque_y'
            o = f'{j}_ground_torque_z'
            labels.append(a)
            labels.append(b)
            labels.append(c)
            labels.append(g)
            labels.append(h)
            labels.append(o)
            labels.append(d)
            labels.append(e)
            labels.append(f)

        with open(f'{i}.mot', 'w') as f:
            f.write('created_by:Wu_chang_fa\n')
            f.write('version=1\n')
            f.write(f'nRows={len(time)}\n')
            f.write(f'nColumns={num_force * 9 + 1}\n')
            f.write('inDegrees=yes\n')
            f.write('endheader\n')
            f.write("time\t" + "\t".join([label for label in labels]))
            f.write("\n")

            for frame_idx in range(frame):
                f.write(f'{time[frame_idx]:.5f}\t')
                for G in range(num_force):
                    f.write(f'{-force[0,(G+1)*6-6,frame_idx]:.5f}\t')
                    f.write(f'{-force[0,(G+1)*6-4,frame_idx]:.5f}\t')
                    f.write(f'{force[0,(G+1)*6-5,frame_idx]:.5f}\t')
                    f.write(f'{force[0,(G+1)*6-3,frame_idx]/-1000:.5f}\t')
                    f.write(f'{-force[0,(G+1)*6-1,frame_idx]/1000:.5f}\t')
                    f.write(f'{-force[0,(G+1)*6-2,frame_idx]/-1000:.5f}\t')
                    f.write(f'{-((force[0, (G + 1) * 6 - 6, frame_idx]*0.052+force[0,(G+1)*6-2,frame_idx])/force[0,(G+1)*6-4,frame_idx]):.5f}\t')
                    f.write(f'{-force[0, (G + 1) * 6 - 2, frame_idx]*0:.5f}\t')
                    f.write(f'{-((-force[0, (G + 1) * 6 - 4, frame_idx]*0.052+force[0,(G+1)*6-3,frame_idx])/force[0,(G+1)*6-4,frame_idx]):.5f}\t')
                f.write('\n')
            f.close()


    elif force_type == 3:
        num_force = int(channel / 8)
        labels = []
        for j in range(num_force):
            a = f'{j}_ground_force_vx'
            b = f'{j}_ground_force_vy'
            c = f'{j}_ground_force_vz'
            d = f'{j}_ground_force_px'
            e = f'{j}_ground_force_py'
            f = f'{j}_ground_force_pz'
            g = f'{j}_ground_torque_x'
            h = f'{j}_ground_torque_y'
            o = f'{j}_ground_torque_z'
            labels.append(a)
            labels.append(b)
            labels.append(c)
            labels.append(g)
            labels.append(h)
            labels.append(o)
            labels.append(d)
            labels.append(e)
            labels.append(f)

        with open(f'{i}.mot', 'w') as f:
            f.write('created_by:Wu_chang_fa\n')
            f.write('version=1\n')
            f.write(f'nRows={len(time)}\n')
            f.write(f'nColumns={num_force * 9 + 1}\n')
            f.write('inDegrees=yes\n')
            f.write('endheader\n')
            f.write("time\t" + "\t".join([label for label in labels]))
            f.write("\n")

            for frame_idx in range(frame):
                f.write(f'{time[frame_idx]:.5f}\t')
                for G in range(num_force):
                    # vx = -(c3d['data']['analogs'][0, (G + 1) * 8 - 8, :].T + c3d['data']['analogs'][0, (G + 1) * 8 - 7,:].T)
                    # f.write(f'{vx[frame_idx]:.5f}\t')
                    f.write(f'{-(force[0,(G+1)*8-8,frame_idx]+force[0,(G+1)*8-7,frame_idx]):.5f}\t')

                    # vy = -(c3d['data']['analogs'][0, (G + 1) * 8 - 4, :].T + c3d['data']['analogs'][0, (G + 1) * 8 - 3,:].T + c3d['data']['analogs'][0,(G + 1) * 8 - 2, :].T + c3d['data']['analogs'][0, (G + 1) * 8 - 1,:].T)
                    # f.write(f'{vy[frame_idx]:.5f}\t')
                    f.write(f'{-(force[0, (G + 1) * 8 - 4, frame_idx] + force[0, (G + 1) * 8 - 3, frame_idx]+ force[0, (G + 1) * 8 - 2, frame_idx]+ force[0, (G + 1) * 8 - 1, frame_idx]):.5f}\t')

                    # vz = (c3d['data']['analogs'][0, (G + 1) * 8 - 6, :].T + c3d['data']['analogs'][0, (G + 1) * 8 - 5,:].T)
                    # f.write(f'{vz[frame_idx]:.5f}\t')
                    f.write(f'{force[0,(G+1)*8-6,frame_idx]+force[0,(G+1)*8-5,frame_idx]:.5f}\t')

                    # tx = (-(c3d['data']['analogs'][0, (G + 1) * 8 - 4, :].T + c3d['data']['analogs'][0, (G + 1) * 8 - 3,:].T - c3d['data']['analogs'][0,(G + 1) * 8 - 2, :].T -c3d['data']['analogs'][0, (G + 1) * 8 - 1, :].T) * 350) / 1000
                    # f.write(f'{tx[frame_idx]:.5f}\t')
                    f.write(f'{(force[0, (G + 1) * 8 - 4, frame_idx] + force[0, (G + 1) * 8 - 3, frame_idx]- force[0, (G + 1) * 8 - 2, frame_idx]- force[0, (G + 1) * 8 - 1, frame_idx])*(-350/1000):.5f}\t')

                    ty = (-(350 * (-c3d['data']['analogs'][0, (G + 1) * 8 - 8, :].T + c3d['data']['analogs'][0,(G + 1) * 8 - 7, :].T) + 210 * (c3d['data']['analogs'][0, (G + 1) * 8 - 6, :].T - c3d['data']['analogs'][0,(G + 1) * 8 - 5,:].T))) / 1000
                    f.write(f'{ty[frame_idx]:.5f}\t')
                    # f.write(f'{(350*(-force[0, (G + 1) * 8 - 8, frame_idx] + force[0, (G + 1) * 8 - 7, frame_idx])+210*(force[0,(G+1)*8-6,frame_idx]+force[0,(G+1)*8-5,frame_idx]))/-1000:.5f}\t')



                    # tz = ((210 * (-c3d['data']['analogs'][0, (G + 1) * 8 - 4, :].T + c3d['data']['analogs'][0,(G + 1) * 8 - 3, :].T +c3d['data']['analogs'][0, (G + 1) * 8 - 2, :].T - c3d['data']['analogs'][0,(G + 1) * 8 - 1, :].T))) / 1000
                    # f.write(f'{tz[frame_idx]:.5f}\t')
                    f.write(f'{(-force[0, (G + 1) * 8 - 4, frame_idx] + force[0, (G + 1) * 8 - 3, frame_idx] +force[0, (G + 1) * 8 - 2, frame_idx] - force[0, (G + 1) * 8 - 1, frame_idx]) * (210/ 1000):.5f}\t')



                    px = -(-(c3d['data']['analogs'][0, (G + 1) * 8 - 8, :].T + c3d['data']['analogs'][0,(G + 1) * 8 - 7, :].T) * (-0.052) - ((210 * (-c3d['data']['analogs'][0, (G + 1) * 8 - 4, :].T + c3d['data']['analogs'][0,(G + 1) * 8 - 3, :].T +c3d['data']['analogs'][0, (G + 1) * 8 - 2, :].T - c3d['data']['analogs'][0,(G + 1) * 8 - 1,:].T))) / 1000) / -(c3d['data']['analogs'][0, (G + 1) * 8 - 4, :].T + c3d['data']['analogs'][0,(G + 1) * 8 - 3, :].T + c3d['data']['analogs'][0, (G + 1) * 8 - 2,:].T +c3d['data']['analogs'][0, (G + 1) * 8 - 1, :].T)
                    f.write(f'{px[frame_idx]:.5f}\t')

                    py = -((c3d['data']['analogs'][0, (G + 1) * 8 - 6, :].T + c3d['data']['analogs'][0, (G + 1) * 8 - 5,:].T) * (-0.052) + (-(c3d['data']['analogs'][0, (G + 1) * 8 - 4, :].T + c3d['data']['analogs'][0,(G + 1) * 8 - 3, :].T - c3d['data']['analogs'][0, (G + 1) * 8 - 2,:].T -c3d['data']['analogs'][0, (G + 1) * 8 - 1, :].T) * 350) / 1000) / -(c3d['data']['analogs'][0, (G + 1) * 8 - 4, :].T + c3d['data']['analogs'][0,(G + 1) * 8 - 3, :].T + c3d['data']['analogs'][0, (G + 1) * 8 - 2,:].T +c3d['data']['analogs'][0, (G + 1) * 8 - 1, :].T)
                    pz = (c3d['data']['analogs'][0, (G + 1) * 8 - 8, :].T + c3d['data']['analogs'][0, (G + 1) * 8 - 7,:].T) * 0
                    f.write(f'{pz[frame_idx]:.5f}\t')
                    f.write(f'{py[frame_idx]:.5f}\t')
                f.write('\n')
            f.close()

    else:
        print('未知测力板类型，暂不支持.mot文件格式转化')



print('文件处理成功：\n')
print('文件顺序为：\n')
for l in range(len(file_paths)):
    print(f'{l}号：{file_paths[l]}\n')


