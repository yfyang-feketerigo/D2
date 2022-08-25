// 20210406 add new pair style: none, for none pair info data file
// class, store configuration generated by lammps 'data' commend
// 20210325 fix displacemnt algorithm, help by yjruan
// 20211123 add members: lx ly lz, automaticlly computed when obj constructed
//			add func: get_lx(), get_ly(), get_lz()
// 20211129 规范namespace
// 20220825 safely cast double to integer
#include "configuration.h"
namespace Configuration
{
	using std::ofstream;
	using std::ifstream;
	using std::clog;
	using std::cout;
	using std::to_string;
	using std::endl;
	using std::cerr;
	using std::string;
	Configuration::Configuration(std::string config_file, BoxType _boxtype, PairStyle _pairstyle, bool _is_sorted)
	{
#ifdef LOG_ON_SCREEN
		clog << "#LAMMPS data file reader..." << '\n';
		clog << "#Ivan Young@CIAC 20201020" << '\n';
		clog << "Reading configuration data file " << config_file << " ..." << '\n';
#endif // LOG_ON_SCREEN
		//timestep = _time;
		filename = config_file;
		std::string firstline;
		ifstream infile;
		HEAD_INFO_LINE = 0;
		infile.open(config_file);
		flag_particle_sorted = _is_sorted;
		if (infile.is_open())
		{
#ifdef LOG_ON_SCREEN
			clog << "Processing data file head information..." << '\n';
#endif // LOG_ON_SCREEN	
			//获取时间步信息
			getline(infile, firstline); //1st line
			string str_timestep;
			for (size_t i = firstline.size() - 1; i > 0; i--)
			{
				if (' ' != firstline[i])
				{
					str_timestep.insert(str_timestep.begin(), firstline[i]);
				}
				else
				{
					break;
				}
			}
			timestep = stoull(str_timestep); //string to unsigned long long

			//跳过描述信息
			infile.ignore(LINE_SKIP_MAX, '\n'); //2nd line

			infile >> particle_num;
			infile.ignore(LINE_SKIP_MAX, '\n'); //3rd line

			infile >> type_num;
			infile.ignore(LINE_SKIP_MAX, '\n'); //4th line

			infile.ignore(LINE_SKIP_MAX, '\n'); //5th line

			HEAD_INFO_LINE = 5;

			/*
			* 读入x边界
			*/
			infile >> xlo >> xhi;
			infile.ignore(LINE_SKIP_MAX, '\n'); //6th line
			HEAD_INFO_LINE++;

			/*
			* 读入y边界
			*/
			infile >> ylo >> yhi;
			infile.ignore(LINE_SKIP_MAX, '\n'); //7th line
			HEAD_INFO_LINE++;

			/*
			* 读入z边界
			*/
			infile >> zlo >> zhi;
			infile.ignore(LINE_SKIP_MAX, '\n'); //8th line
			HEAD_INFO_LINE++;


			if (_boxtype == BoxType::tilt) //三斜盒子的处理
			{
				infile >> xy >> xz >> yz;
				infile.ignore(LINE_SKIP_MAX, '\n'); //9th line
				HEAD_INFO_LINE++;
			}
			infile.ignore(LINE_SKIP_MAX, '\n'); //10th line
			HEAD_INFO_LINE++;


			// setting lx ly lz
			lx = xhi - xlo;
			ly = yhi - ylo;
			lz = zhi - zlo;


			string string_mass_info; //读入粒子质量信息
			for (size_t i = 0; i < 3 + type_num; i++)
			{
				getline(infile, string_mass_info);
				strvec_mass_info.push_back(string_mass_info);
			} //11th - 14th line
			HEAD_INFO_LINE += (3 + type_num);
			//infile.ignore(LINE_SKIP_MAX, '\n'); //15th line

			size_t total_pair_line = 0;
			size_t pair_info_space_line = 3;

			switch (_pairstyle) //处理粒子势能信息
			{
			case Configuration::PairStyle::single:
				total_pair_line = type_num;
				break;
			case Configuration::PairStyle::pair:
				total_pair_line += type_num;
				total_pair_line += type_num * (type_num - 1) / 2;
				break;
			case Configuration::PairStyle::none:
				total_pair_line = 0;
				pair_info_space_line = 0;
			default:
				break;
			}
			string string_pair_info;
			for (size_t i = 0; i < pair_info_space_line + total_pair_line; i++)
			{
				getline(infile, string_pair_info);
				strvec_pair_info.push_back(string_pair_info);
			} //15th - 19th line
			HEAD_INFO_LINE += strvec_pair_info.size();


			getline(infile, str_atoms_info); //粒子描述信息
			HEAD_INFO_LINE++;
			infile.ignore(LINE_SKIP_MAX, '\n');
			HEAD_INFO_LINE++;
#ifdef LOG_ON_SCREEN
			clog << "Head information has been processed" << '\n';
#endif // LOG_ON_SCREEN
}
		else
		{
			cerr << "File " << config_file << " open failed!" << endl;
			throw std::runtime_error(config_file.c_str());
		}
		infile.close();

		/*
		* 屏幕输出读入的data文件信息，用于校验差错
		*/
#ifdef LOG_ON_SCREEN
		clog << firstline << '\n';;
		clog << "Configuration data file " << config_file << " has " << particle_num << " particles" << '\n';
		clog << "Configuration has " << type_num << " particle type(s)" << '\n';
		clog << "Time Step: " << timestep << '\n';
		for (size_t i = 0; i < strvec_mass_info.size(); i++)
		{
			clog << strvec_mass_info[i] << '\n';
		}

		for (size_t i = 0; i < strvec_pair_info.size(); i++)
		{
			clog << strvec_pair_info[i] << '\n';
		}
		clog << '\n' << "Box Parameters:" << '\n';
		clog << "xlo, xhi: " << xlo << ' ' << xhi << '\n';
		clog << "ylo, yhi: " << ylo << ' ' << yhi << '\n';
		clog << "zlo, zhi: " << zlo << ' ' << zhi << '\n';
		clog << "xy, xz, yz: " << xy << ' ' << xz << ' ' << yz << '\n';
		//clog << "head line info read!" << endl;
		clog << '\n';
		clog << str_atoms_info << '\n';
		clog << '\n' << "File HEAD LINE: " << HEAD_INFO_LINE << '\n';
		clog << "File GAP LINE: " << GAP_LINE << '\n';
		clog << '\n';
		clog << "Reading coordinates..." << '\n';
#endif // LOG_ON_SCREEN
		Input in_data(config_file, HEAD_INFO_LINE);
		/*
		* 开始处理粒子坐标信息
		*/
		in_data.open_file();
		in_data.skiphead();
		vec_particle.resize(particle_num);
		for (size_t i = 0; i < particle_num; i++)
		{
			in_data.read_line_data();
			vec_particle[i].id = std::lround(in_data.get_data()[0]);
			vec_particle[i].type = std::rint(in_data.get_data()[1]);
			vec_particle[i].rx = in_data.get_data()[2];
			vec_particle[i].ry = in_data.get_data()[3];
			vec_particle[i].rz = in_data.get_data()[4];
			vec_particle[i].box_x = std::rint(in_data.get_data()[5]);
			vec_particle[i].box_y = std::rint(in_data.get_data()[6]);
			vec_particle[i].box_z = std::rint(in_data.get_data()[7]);
		}

		in_data.skip_line(GAP_LINE); //跳过坐标与速度间空行
#ifdef LOG_ON_SCREEN
		clog << "Coordinates have been read!" << '\n';
		clog << "Reading velocities..." << '\n';
#endif // LOG_ON_SCREEN
		/*
		* 开始处理粒子速度信息
		*/
		for (size_t i = 0; i < particle_num; i++)
		{
			in_data.read_line_data();
			Particle& p_particle = seek_id(vec_particle, (size_t)in_data.get_data()[0]);
			p_particle.vx = in_data.get_data()[1];
			p_particle.vy = in_data.get_data()[2];
			p_particle.vz = in_data.get_data()[3];
		}
#ifdef LOG_ON_SCREEN
		clog << "Velocities have been read!" << '\n';
		clog << "Configuration data file " << config_file << " has been read!" << '\n';
		clog << '\n';
#endif // LOG_ON_SCREEN
		infile.close();
		}

	Configuration::Configuration(const Configuration& conf)
	{
		this->HEAD_INFO_LINE = conf.HEAD_INFO_LINE;
		this->particle_num = conf.particle_num;
		this->timestep = conf.timestep;
		this->xlo = conf.xlo;
		this->ylo = conf.ylo;
		this->zlo = conf.zlo;
		this->xhi = conf.xhi;
		this->yhi = conf.yhi;
		this->zhi = conf.zhi;
		this->xy = conf.xy;
		this->xz = conf.xz;
		this->yz = conf.yz;
		this->lx = conf.lx;
		this->ly = conf.ly;
		this->lz = conf.lz;
		this->type_num = conf.type_num;
		this->strvec_mass_info = conf.strvec_mass_info;
		this->strvec_pair_info = conf.strvec_pair_info;
		this->str_atoms_info = conf.str_atoms_info;
		this->vec_particle = conf.vec_particle;
		if (conf.flag_particle_sorted)
		{
			this->sort_particle();
		}
		else
		{
			this->flag_particle_sorted = false;
			this->pvec_particle_sorted.clear();
			this->pvec_particle_sorted = std::move(std::vector<const Particle*>());
		}
	}

	Configuration& Configuration::operator=(const Configuration& conf)
	{
		if (this != &conf)
		{
			this->HEAD_INFO_LINE = conf.HEAD_INFO_LINE;
			this->particle_num = conf.particle_num;
			this->timestep = conf.timestep;
			this->xlo = conf.xlo;
			this->ylo = conf.ylo;
			this->zlo = conf.zlo;
			this->xhi = conf.xhi;
			this->yhi = conf.yhi;
			this->zhi = conf.zhi;
			this->xy = conf.xy;
			this->xz = conf.xz;
			this->yz = conf.yz;
			this->lx = conf.lx;
			this->ly = conf.ly;
			this->lz = conf.lz;
			this->type_num = conf.type_num;
			this->strvec_mass_info = conf.strvec_mass_info;
			this->strvec_pair_info = conf.strvec_pair_info;
			this->str_atoms_info = conf.str_atoms_info;
			this->vec_particle = conf.vec_particle;
			if (conf.flag_particle_sorted)
			{
				this->sort_particle();
			}
			else
			{
				this->flag_particle_sorted = false;
				this->pvec_particle_sorted.clear();
				this->pvec_particle_sorted = std::move(std::vector<const Particle*>());
			}
		}
		return *this;
	}

	const Particle& Configuration::get_particle(size_t _id) const //return particle with given ID
	{
		for (size_t i = 0; i < vec_particle.size(); i++)
		{
			if (_id == vec_particle[i].id)
			{
				return vec_particle[i];
			}
		}
		std::cerr << "particle " << _id << " not found";
		throw std::runtime_error("particle " + std::to_string(_id) + " not found");
	}

	void Configuration::sort_particle()
	{
		auto comp_id = [](const Particle* a, const Particle* b) {return a->id <= b->id; };
		pvec_particle_sorted.resize(vec_particle.size());
		for (size_t i = 0; i < pvec_particle_sorted.size(); i++)
		{
			pvec_particle_sorted[i] = &vec_particle[i];
		}
		std::sort(pvec_particle_sorted.begin(), pvec_particle_sorted.end(), comp_id);
		flag_particle_sorted = true;
	}

	const std::vector<const Particle*>& Configuration::get_pvec_particle_sorted() const
	{
#ifdef SAFE_GET_PARTICLE_SORTED
		if (!flag_particle_sorted)
			throw std::runtime_error("particle not sorted when using get_particle_sorted()!!!");
#endif // SAFE_GET_PARTICLE_SORTED
		return pvec_particle_sorted;
	}

	const Particle& Configuration::get_particle_sorted(size_t id) const // CAUTION! this method need all particle id be consecutive and start with 0;
	{
		//std::cout << "using sorted!" << '\n';
#ifdef SAFE_GET_PARTICLE_SORTED
		if (!flag_particle_sorted)
			throw std::runtime_error("particle not sorted when using get_particle_sorted(size_t id)!!!");
		if (id > get_particle().size() || id <= 0)
			throw std::runtime_error(("Searching particle id " + std::to_string(id) + " illegal\n").c_str());
		if (id != pvec_particle_sorted[id - 1]->id)
			throw std::runtime_error(("id not match when searching particle sorted! this could due to inconsecutive id.\n"));
#endif // SAFE_GET_PARTICLE_SORTED
		return *(pvec_particle_sorted[id - 1]);
	}

	size_t Configuration::__add_particle(const Particle& new_pa) //新添加一个粒子
	{
		bool flag_inbox = new_pa.rx >= xlo && new_pa.rx <= xhi
			&& new_pa.ry >= ylo && new_pa.ry <= yhi
			&& new_pa.rz >= zlo && new_pa.rz <= zhi;
		if (!flag_inbox)
		{
			cerr << "new particle coordiantion is not in box!" << endl;
			throw std::runtime_error("new particle coordiantion is not in box!");
		}
		bool flag_oldtype = true;
		for (size_t i = 0; i < vec_particle.size(); i++)
		{
			flag_oldtype = flag_oldtype && (new_pa.type == vec_particle[i].type);
		}
		if (!flag_oldtype)
		{
			cerr << "WARNING: NEW particle type add, make sure pair&mass info be modified" << endl;
		}
		particle_num++;
		vec_particle.push_back(new_pa);
		flag_particle_sorted = false;
		return vec_particle.size();
	}

	void Configuration::to_data(string fname, BoxType _boxtype) //以lammps data文件格式输出
	{
		ofstream ofile;
		ofile.open(fname);
		if (!ofile.is_open())
		{
			cerr << fname << " open failed" << endl;
			throw std::runtime_error((fname + " open failed").c_str());
		}
		ofile << "LAMMPS data file via C++, Configuration class, timestep = " << timestep << '\n';
		ofile << '\n';
		ofile << particle_num << " atoms" << '\n';
		ofile << type_num << " atom types" << '\n';
		ofile << '\n';

		ofile << xlo << " " << xhi << " " << "xlo " << "xhi" << '\n';
		ofile << ylo << " " << yhi << " " << "ylo " << "yhi" << '\n';
		ofile << zlo << " " << zhi << " " << "zlo " << "zhi" << '\n';
		if (_boxtype == BoxType::tilt)
		{
			ofile << xy << " " << xz << " " << yz << " " << "xy " << "xz " << "yz";
		}
		for (size_t i = 0; i < strvec_mass_info.size(); i++)
		{
			ofile << strvec_mass_info[i] << '\n';
		}
		for (size_t i = 0; i < strvec_pair_info.size(); i++)
		{
			ofile << strvec_pair_info[i] << '\n';
		}
		ofile << str_atoms_info << '\n';
		ofile << '\n';
		for (size_t i = 0; i < vec_particle.size(); i++)
		{
			Particle& pa = vec_particle[i];
			ofile << pa.id << " " << pa.type << " " << pa.rx << " " << pa.ry << " " << pa.rz << " "
				<< pa.box_x << " " << pa.box_y << " " << pa.box_z << '\n';
		}
		ofile << '\n';
		ofile << "Velocities" << '\n';
		ofile << '\n';
		for (size_t i = 0; i < vec_particle.size(); i++)
		{
			Particle& pa = vec_particle[i];
			ofile << pa.id << " " << pa.vx << " " << pa.vy << " " << pa.vz << '\n';
		}
		ofile.close();
		return;
	}

	void Configuration::to_dump(string ofname, string opath, string style) const //以lammps dumpo文件格式输出, style 可以指定为不同样式
	{
#ifdef LOG_ON_SCREEN
		clog << "converting data file to dump file..." << '\n';
		clog << "dump style: " << style << '\n';
		clog << "dump path: " << opath << '\n';
		clog << "dump file name: " << ofname << '\n';
#endif // LOG_ON_SCREEN
		if (!boost::filesystem::exists(opath))
		{
			boost::filesystem::create_directories(opath);
		}
		string full_opath = opath + ofname;

		ofstream ofile;
		ofile.open(full_opath);
		if (!ofile.is_open())
		{
			//cerr << full_opath << " open failed" << endl;
			throw std::runtime_error((full_opath + " open failed").c_str());
		}
		/*
		double lx = xhi - xlo;
		double ly = yhi - ylo;
		double lz = zhi - zlo;
		*/
		if (style == "yjruan")
		{
			ofile << "ITEM: TIMESTEP" << '\n';
			ofile << timestep << '\n';
			ofile << "ITEM: NUMBER OF ATOMS" << '\n';
			ofile << particle_num << '\n';
			ofile << "ITME: BOX BOUNDS xy xz yz pp pp pp" << '\n';
			ofile.precision(15);

			/*
			* 注意此处边界的处理，data文件与dump文件格式不相同
			*/
			double xlo_bound = xlo + std::min({ 0.0, xy, xz, xy + xz });
			double xhi_bound = xhi + std::max({ 0.0, xy, xz, xy + xz });

			double ylo_bound = ylo + std::min({ 0.0, yz });
			double yhi_bound = yhi + std::max({ 0.0, yz });

			double zlo_bound = zlo;
			double zhi_bound = zhi;

			ofile << xlo_bound << " " << xhi_bound << " " << xy << '\n';
			ofile << ylo_bound << " " << yhi_bound << " " << xz << '\n';
			ofile << zlo_bound << " " << zhi_bound << " " << yz << '\n';
			ofile << "ITEM: ATOMES id xu yu zu" << '\n';
			for (size_t i = 0; i < vec_particle.size(); i++)
			{
				const Particle& p_pa = vec_particle[i];
				double xu = get_rx_real(p_pa);
				double yu = get_ry_real(p_pa);
				double zu = get_rz_real(p_pa);
				size_t id = p_pa.id;
				ofile << id << " ";
				ofile << xu << " ";
				ofile << yu << " ";
				ofile << zu << " ";
				ofile << '\n';
			}
			return;
		}

		cerr << "wrong dump style: " << style << '\n';
		throw std::runtime_error(("wrong dump style: " + style).c_str());

	}
}