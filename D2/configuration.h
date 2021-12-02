// 20210406 add new pair style: none, for none pair info data file
// class, store configuration generated by lammps 'data' commend
// 20210325 fix displacemnt algorithm, help by yjruan
// 20211123 add members: lx ly lz, automaticlly computed when obj constructed
//			add func: get_lx(), get_ly(), get_lz()
// 20211129 规范namespace
/*
* 用于读取KA模型（纯LJ粒子）体系LAMMPS data文件
* 同时允许输出单粒子的参量，使用para_to_dump方法
*/
#pragma once
#include "particle.h"
#include "input.h"
#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <array>
#include <vector>
#include <initializer_list>
#include <algorithm>
//#define LOG_ON_SCREEN
namespace Configuration
{

	class Configuration
	{
	private:
		static const auto LINE_SKIP_MAX = std::numeric_limits<std::streamsize>::max();
		//static const size_t LINE_SKIP_MAX = 2147483647;
		static const size_t GAP_LINE = 3;		 //description line between position and velocity

		size_t HEAD_INFO_LINE = 0; //data file description information, e.g. mass, pair etc.
		size_t particle_num = 0; //total particle number
		unsigned long long timestep = 0;

		double xlo = 0;//xlo difined in lammps
		double ylo = 0;//ylo defined in lammps
		double zlo = 0;//zlo defined in lammps

		double xhi = 0;//xhi difined in lammps
		double yhi = 0;//yhi defined in lammps
		double zhi = 0;//zhi defined in lammps

		double xy = 0;//xy defined in lammps
		double xz = 0;//xz defined in lammps
		double yz = 0;//yz defined in lammps

		double lx = 0;
		double ly = 0;
		double lz = 0;

		size_t type_num = 0;//number of particle types

		std::vector<std::string> strvec_mass_info; //store mass information in data file
		std::vector<std::string> strvec_pair_info; //store pair information in data file
		std::string str_atoms_info;
		std::vector<Particle> vec_particle; //particles container
		std::string filename;

	public:
		enum class BoxType
		{
			orthogonal,
			tilt
		};
		enum class PairStyle
		{
			single,
			pair,
			none
		};
		Configuration(std::string config_file, BoxType _boxtype = BoxType::orthogonal, PairStyle _pairstyle = PairStyle::single);				   //
		Configuration() {};
		inline size_t GET_LINE_MAX() { return LINE_SKIP_MAX; }
		inline size_t GET_HEAD_INFO_LINE() { return HEAD_INFO_LINE; }
		inline static size_t GET_GAP_LINE() { return GAP_LINE; }


		inline void SET_HEAD_INFO_LINE(size_t _HEAD_INFO_LINE) { HEAD_INFO_LINE = _HEAD_INFO_LINE; }
		//inline static void SET_GAP_LINE(size_t _GAP_LINE) { GAP_LINE = _GAP_LINE; }

		inline void set_time_step(unsigned long long _timestep) { timestep = _timestep; };

		inline void set_particle_num(size_t _num) { particle_num = _num; } //set total number of particles
		inline size_t get_particle_num() const { return particle_num; }	   //return total number of particles
		inline size_t get_type_num() { return type_num; }				   //return number of particle types
		inline const std::vector<std::string>& get_mass_info() const				//return mass info
		{
			return strvec_mass_info;
		}
		inline const std::vector<std::string>& get_pair_info() const			//return pair info
		{
			return strvec_pair_info;
		}

		inline double get_xlo() const { return xlo; }//return xhi defined in lammps
		inline double get_xhi() const { return xhi; }//return xhi defined in lammps
		inline double get_ylo() const { return ylo; }//return ylo defined in lammps
		inline double get_yhi() const { return yhi; }//return yhi defined in lammps
		inline double get_zlo() const { return zlo; }//return zlo defined in lammps
		inline double get_zhi() const { return zhi; }//return zhi defined in lammps
		inline double get_xy() const { return xy; }//return xy defined in lammps
		inline double get_yz() const { return yz; }//return yz defined in lammps
		inline double get_xz() const { return xz; }//return xz defined in lammps
		inline double get_lx() const { return lx; }
		inline double get_ly() const { return ly; }
		inline double get_lz() const { return lz; }

		/*
		* ri_real部分用于获得粒子真实坐标，以方便计算位移
		*/
		inline double get_rx_real(const Particle& pa) const
		{
			double rx_real = pa.rx + pa.box_x * (xhi - xlo)
				+ pa.box_y * xy + pa.box_z * xz;
			return rx_real;
		}

		inline double get_ry_real(const Particle& pa) const
		{
			double ry_real = pa.ry + pa.box_y * (yhi - ylo)
				+ pa.box_z * yz;
			return ry_real;
		}

		inline double get_rz_real(const Particle& pa) const
		{
			double rz_real = pa.rz + pa.box_z * (zhi - zlo);
			return rz_real;
		}

		/*
		* 手动修改边界
		*/

		inline double set_xboundary(double _xlo, double _xhi)
		{
			xlo = _xlo;
			xhi = _xhi;
		}
		inline double set_yboundary(double _ylo, double _yhi)
		{
			ylo = _ylo;
			yhi = _yhi;
		}
		inline double set_zboundary(double _zlo, double _zhi)
		{
			zlo = _zlo;
			zhi = _zhi;
		}
		inline double set_tilt(double _xy, double _xz, double _yz)
		{
			xy = _xy;
			xz = _xz;
			yz = _yz;
		}

		inline unsigned long long get_timestep() const { return timestep; }

		inline const std::string& get_filename() const { return filename; }

		inline const Particle& get_particle(size_t _id) const //return particle with given ID
		{
			for (size_t i = 0; i < vec_particle.size(); i++)
			{
				if (_id == vec_particle[i].id)
				{
					return vec_particle[i];
				}
			}
			std::cerr << "particle " << _id << " not found";
			throw("particle " + std::to_string(_id) + " not found");
		}
		inline const std::vector<Particle>& get_particle() const//return vector of all particles
		{
			return vec_particle;
		}

		size_t __add_particle(const Particle& new_pa);
		inline void __clear_vec_particle()
		{
			vec_particle.clear();
			particle_num = 0;
		}
		void to_data(std::string fname, BoxType _boxtype = BoxType::tilt);

		template<typename T>
		void para_to_dump(std::string fname, std::initializer_list<std::string> add_para_name, std::initializer_list<std::vector<T>> add_para, std::vector<std::string> comments = {}) const//注意额外参量与粒子序号的对应关系
		{
			std::ofstream ofile;
			ofile.open(fname);
			if (!ofile.is_open())
			{
				std::cerr << fname << " open failed" << std::endl;
				throw (fname + " open failed");
			}
			for (size_t i = 0; i < comments.size(); i++)
			{
				ofile << comments[i] << '\n';
			}
			ofile << "ITEM: TIMESTEP" << '\n';
			ofile << timestep << '\n';
			ofile << "ITEM: NUMBER OF ATOMS" << '\n';
			ofile << particle_num << '\n';
			ofile << "ITEM: BOX BOUNDS xy xz yz pp pp pp " << '\n';
			/*
			* 注意这一部分输出的体系边界方式是按照ovito识别的方式
			* 这一输出方式与data文件不同，但与dump指令相同
			*/
			auto xvi = { 0.,xy,xz,xy + xz };
			auto x_minmax = std::minmax_element(xvi.begin(), xvi.end());
			double visual_xlo = xlo + *x_minmax.first;
			double visual_xhi = xhi + *x_minmax.second;
			auto yvi = { 0.,yz };
			auto y_minmax = std::minmax_element(yvi.begin(), yvi.end());
			double visual_ylo = ylo + *y_minmax.first;
			double visual_yhi = yhi + *y_minmax.second;
			ofile << visual_xlo << " " << visual_xhi << " " << xy << '\n';
			ofile << visual_ylo << " " << visual_yhi << " " << xz << '\n';
			ofile << zlo << " " << zhi << " " << yz << '\n';
			ofile << "ITEM: ATOMS id type x y z ix iy iz ";
			//输出额外参量名
			for (auto it_li = add_para_name.begin(); it_li < add_para_name.end(); it_li++)
			{
				ofile << *it_li << " ";
			}
			ofile << '\n';
			//输出粒子信息，随后输出额外参量值
			for (size_t i = 0; i < vec_particle.size(); i++)
			{
				const Particle& pa = vec_particle[i];
				ofile << pa.id << " " << pa.type << " " << pa.rx << " " << pa.ry << " " << pa.rz << " ";
				ofile << pa.box_x << " " << pa.box_y << " " << pa.box_z << " ";
				for (auto it_li = add_para.begin(); it_li < add_para.end(); it_li++)
				{
					ofile << (*it_li)[i] << " ";
				}
				ofile << '\n';
			}
			ofile.close();
			return;
		}

		void to_dump(std::string ofname, std::string opath = "./", std::string style = "") const;
	};
}
