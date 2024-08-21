#include "lammps_dump.h"
#include "string_tools.h"
#include <fmt/format.h>
#include <iostream>


LMP_dump::LMP_dump(std::string path) {
    std::fstream fid;
    fid.open(path, std::fstream::in);
    if (!fid) {
        std::cout << "Dump file does not exist.\n";
        exit(0);
    }
    while(fid) {
        auto line = read_line(fid);
        auto args = split(line);
        if (args.empty()) continue;
        else if (line == "ITEM: TIMESTEP") {
            int n = from_string<int>(read_line(fid));
            timesteps.push_back(n);
        }
        else if (line == "ITEM: NUMBER OF ATOMS") {
            num_atoms = from_string<int>(read_line(fid));
        }
        else if (startswith(line, "ITEM: BOX BOUNDS")) {
            _timestep_box = MatrixXd::Zero(3, 3);
            double xlo, xhi, ylo, yhi, zlo, zhi;
            for (auto arg : args) {
                if (arg == "xy" || arg == "xz" || arg == "yz") {
                    triclinic = true;
                }
            }
            if (triclinic == true) {
                // Triclinic box : https://lammps.sandia.gov/doc/Howto_triclinic.html
                args = split(read_line(fid));
                auto xlo_bound = str2dbl(args[0]);
                auto xhi_bound = str2dbl(args[1]);
                auto xy = str2dbl(args[2]);
                args = split(read_line(fid));
                auto ylo_bound = str2dbl(args[0]);
                auto yhi_bound = str2dbl(args[1]);
                auto xz = str2dbl(args[2]);
                args = split(read_line(fid));
                auto zlo_bound = str2dbl(args[0]);
                auto zhi_bound = str2dbl(args[1]);
                auto yz = str2dbl(args[2]);

                xlo = xlo_bound - std::min(0.0, std::min(xy, std::min(xz, xy+xz)));
                xhi = xhi_bound - std::max(0.0, std::max(xy, std::max(xz, xy+xz)));
                ylo = ylo_bound - std::min(0.0, yz);
                yhi = yhi_bound - std::max(0.0, yz);
                zlo = zlo_bound;
                zhi = zhi_bound;

                _timestep_box(0, 0) = xhi - xlo;
                _timestep_box(0, 1) = xy;
                _timestep_box(0, 2) = xz;
                _timestep_box(1, 0) = 0.0;
                _timestep_box(1, 1) = yhi - ylo;
                _timestep_box(1, 2) = yz;
                _timestep_box(2, 0) = 0.0;
                _timestep_box(2, 1) = 0.0;
                _timestep_box(2, 2) = zhi - zlo;
            }
            else {
                args = split(read_line(fid));
                xlo = str2dbl(args[0]);
                xhi = str2dbl(args[1]);
                args = split(read_line(fid));
                ylo = str2dbl(args[0]);
                yhi = str2dbl(args[1]);
                args = split(read_line(fid));
                zlo = str2dbl(args[0]);
                zhi = str2dbl(args[1]);

                _timestep_box(0, 0) = xhi - xlo;
                _timestep_box(1, 1) = yhi - ylo;
                _timestep_box(2, 2) = zhi - zlo;
            }
            box.push_back(_timestep_box);
        }
        else if (startswith(line, "ITEM: ATOMS")) {
            if (std::find(args.begin(), args.end(), "xs") != args.end()) {
                coords_type = "scaled";
            }
            else if (std::find(args.begin(), args.end(), "x") != args.end()) {
                coords_type = "unscaled";
            }
            else if (std::find(args.begin(), args.end(), "xu") != args.end()) {
                coords_type = "unwrapped";
            }
            else {
                fmt::print("Unrecognized coordinate type.\n");
                fmt::print("Recognized coordinate types are x, xs, and xu.\n");
                exit(0);
            }
            MatrixXd timestep_coords(num_atoms, 3);
            int counter = 0;
            while(fid && counter < num_atoms) {
                line = read_line(fid);
                args = split(line);
                int num_args = args.size();
                int row_no = from_string<int>(args[0]) - 1;
                timestep_coords(row_no, 0) = str2dbl(args[num_args - 3]);
                timestep_coords(row_no, 1) = str2dbl(args[num_args - 2]);
                timestep_coords(row_no, 2) = str2dbl(args[num_args - 1]);
                counter++;
            }
            atom_coords.push_back(timestep_coords);
        }
    }
    fmt::print("Read {} timesteps from {}.\n", timesteps.size(), path);
}

