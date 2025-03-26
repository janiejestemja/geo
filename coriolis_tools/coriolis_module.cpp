#include <Python.h>
#include <cmath>
#include <vector>


// Function to calculate Earth's radius at a given latitude
double radius_at_latitude(double lat) {
    double lat_rad = lat * M_PI / 180.0; // Convert latitude from degrees to radians
    return 6371.0 * (1 - 0.5 * (1 - std::cos(lat_rad))) + 
           6378.137 * (0.5 * (1 - std::cos(lat_rad)));
}

// Wrapper for radius_at_latitude
static PyObject* py_radius_at_latitude(PyObject* self, PyObject* args) {
    double lat;
    if (!PyArg_ParseTuple(args, "d", &lat)) {
        return nullptr;
    }
    double radius = radius_at_latitude(lat);
    return Py_BuildValue("d", radius);
}

// Function to calculate the great-circle distance using the Haversine formula
double haversine(double lat1, double lon1, double lat2, double lon2) {
    // Convert latitude and longitude from degrees to radians
    lat1 = lat1 * M_PI / 180.0;
    lon1 = lon1 * M_PI / 180.0;
    lat2 = lat2 * M_PI / 180.0;
    lon2 = lon2 * M_PI / 180.0;

    // Haversine formula
    double dlon = lon2 - lon1;
    double dlat = lat2 - lat1;
    double a = std::pow(std::sin(dlat / 2.0), 2) + 
               std::cos(lat1) * std::cos(lat2) * std::pow(std::sin(dlon / 2.0), 2);
    double c = 2.0 * std::atan2(std::sqrt(a), std::sqrt(1.0 - a));

    // Get Earth's radius based on the average latitude of the two points
    double avg_lat = (lat1 + lat2) / 2.0 * 180.0 / M_PI;  // Convert back to degrees for radius calc
    double r = radius_at_latitude(avg_lat);

    // Distance in kilometers
    return c * r;
}

// Wrapper for the haversine function
static PyObject* py_haversine(PyObject* self, PyObject* args) {
    double lat1, lon1, lat2, lon2;
    if (!PyArg_ParseTuple(args, "dddd", &lat1, &lon1, &lat2, &lon2)) {
        return nullptr;
    }
    double distance = haversine(lat1, lon1, lat2, lon2);
    return Py_BuildValue("d", distance);
}


// Function to calculate the direction vector between two lat/lon points
std::vector<double> direction_vector(double lat1, double lon1, double lat2, double lon2) {
    lat1 = lat1 * M_PI / 180.0;
    lon1 = lon1 * M_PI / 180.0;
    lat2 = lat2 * M_PI / 180.0;
    lon2 = lon2 * M_PI / 180.0;

    double x1 = std::cos(lat1) * std::cos(lon1);
    double y1 = std::cos(lat1) * std::sin(lon1);
    double z1 = std::sin(lat1);

    double x2 = std::cos(lat2) * std::cos(lon2);
    double y2 = std::cos(lat2) * std::sin(lon2);
    double z2 = std::sin(lat2);

    std::vector<double> direction = {x2 - x1, y2 - y1, z2 - z1};

    double norm = std::sqrt(direction[0] * direction[0] + direction[1] * direction[1] + direction[2] * direction[2]);
    direction[0] /= norm;
    direction[1] /= norm;
    direction[2] /= norm;

    return direction;
}

// Wrapper for direction_vector
static PyObject* py_direction_vector(PyObject* self, PyObject* args) {
    double lat1, lon1, lat2, lon2;
    if (!PyArg_ParseTuple(args, "dddd", &lat1, &lon1, &lat2, &lon2)) {
        return nullptr;
    }
    std::vector<double> direction = direction_vector(lat1, lon1, lat2, lon2);
    return Py_BuildValue("ddd", direction[0], direction[1], direction[2]);
}

// Function to calculate the rotation matrix based on latitude
std::vector<std::vector<double>> rotation_matrix(double latitude) {
    double lat_rad = latitude * M_PI / 180.0; // Convert latitude to radians
    const double omega = 7.2921159e-5; // rad/s

    return {
        {0, -omega * std::cos(lat_rad), omega * std::sin(lat_rad)},
        {omega * std::cos(lat_rad), 0, 0},
        {-omega * std::sin(lat_rad), 0, 0}
    };
}

// Wrapper for rotation_matrix
static PyObject* py_rotation_matrix(PyObject* self, PyObject* args) {
    double latitude;
    if (!PyArg_ParseTuple(args, "d", &latitude)) {
        return nullptr;
    }
    auto matrix = rotation_matrix(latitude);
    return Py_BuildValue("ddd(ddd)", 
        matrix[0][0], matrix[0][1], matrix[0][2],
        matrix[1][0], matrix[1][1], matrix[1][2],
        matrix[2][0], matrix[2][1], matrix[2][2]);
}

// Coriolis acceleration calculation
std::vector<std::vector<double>> coriolis_acc(double lat1, double lat2, double x_dir, double y_dir, double z_dir, double time, int num_steps) {
    std::vector<std::vector<double>> coriolis_accelerations(num_steps, std::vector<double>(3));
    
    for (int i = 0; i < num_steps; ++i) {
        double t = (time * i) / (num_steps - 1);
        double current_latitude = lat1 + (t / time) * (lat2 - lat1);
        
        auto R = rotation_matrix(current_latitude);
        // Calculate Coriolis acceleration
        for (int j = 0; j < 3; ++j) {
            coriolis_accelerations[i][j] = R[j][0] * x_dir +
                                           R[j][1] * y_dir +
                                           R[j][2] * z_dir;
        }
    }
    return coriolis_accelerations;
}

// Python wrapper for coriolis_acc
static PyObject* py_coriolis_acc(PyObject* self, PyObject* args) {
    double lat1, lat2, x_dir, y_dir, z_dir, time;
    int num_steps;

    // Parse Python arguments: lat1, lat2, x_dir, y_dir, z_dir, time, num_steps
    if (!PyArg_ParseTuple(args, "ddddddi", &lat1, &lat2, &x_dir, &y_dir, &z_dir, &time, &num_steps)) {
        return nullptr;
    }

    // Call the C++ coriolis_acc function
    std::vector<std::vector<double>> result = coriolis_acc(lat1, lat2, x_dir, y_dir, z_dir, time, num_steps);

    // Convert the result to a Python list of lists
    PyObject* py_result = PyList_New(num_steps);
    for (int i = 0; i < num_steps; ++i) {
        PyObject* py_row = PyList_New(3);
        PyList_SetItem(py_row, 0, PyFloat_FromDouble(result[i][0]));
        PyList_SetItem(py_row, 1, PyFloat_FromDouble(result[i][1]));
        PyList_SetItem(py_row, 2, PyFloat_FromDouble(result[i][2]));
        PyList_SetItem(py_result, i, py_row);
    }

    return py_result;
}


// Calculate velocities based on Coriolis acceleration
std::vector<std::vector<double>> calculate_velocities(const std::vector<std::vector<double>>& coriolis_acceleration, double time, int num_steps) {
    std::vector<std::vector<double>> coriolis_velocity(num_steps, std::vector<double>(3, 0.0));

    for (int i = 1; i < num_steps; ++i) {
        double dt = time / (num_steps - 1);
        for (int j = 0; j < 3; ++j) {
            coriolis_velocity[i][j] = coriolis_velocity[i - 1][j] + 
                0.5 * (coriolis_acceleration[i][j] + coriolis_acceleration[i - 1][j]) * dt;
        }
    }
    return coriolis_velocity;
}

// Wrapper for calculate_velocities
static PyObject* py_calculate_velocities(PyObject* self, PyObject* args) {
    PyObject* py_coriolis_acceleration;
    double time;
    int num_steps;

    if (!PyArg_ParseTuple(args, "OdI", &py_coriolis_acceleration, &time, &num_steps)) {
        return nullptr;
    }

    // Extract coriolis_acceleration from Python list of lists
    std::vector<std::vector<double>> coriolis_acceleration;
    if (PyList_Check(py_coriolis_acceleration)) {
        for (Py_ssize_t i = 0; i < PyList_Size(py_coriolis_acceleration); ++i) {
            PyObject* inner_list = PyList_GetItem(py_coriolis_acceleration, i);
            std::vector<double> inner_vector;

            if (PyList_Check(inner_list)) {
                for (Py_ssize_t j = 0; j < PyList_Size(inner_list); ++j) {
                    inner_vector.push_back(PyFloat_AsDouble(PyList_GetItem(inner_list, j)));
                }
            }
            coriolis_acceleration.push_back(inner_vector);
        }
    }

    auto velocities = calculate_velocities(coriolis_acceleration, time, num_steps);

    // Build Python list of lists for output
    PyObject* py_result = PyList_New(num_steps);
    for (size_t i = 0; i < velocities.size(); ++i) {
        PyObject* py_inner = PyList_New(3);
        for (size_t j = 0; j < 3; ++j) {
            PyList_SetItem(py_inner, j, Py_BuildValue("d", velocities[i][j]));
        }
        PyList_SetItem(py_result, i, py_inner);
    }
    
    return py_result;
}

// Calculate drift distances based on velocities
std::vector<std::vector<double>> calculate_drift_distances(const std::vector<std::vector<double>>& coriolis_velocity, double time, int num_steps) {
    std::vector<std::vector<double>> coriolis_drift_distance(num_steps, std::vector<double>(3, 0.0));

    for (int i = 1; i < num_steps; ++i) {
        double dt = time / (num_steps - 1);
        for (int j = 0; j < 3; ++j) {
            coriolis_drift_distance[i][j] = coriolis_drift_distance[i - 1][j] + 
                0.5 * (coriolis_velocity[i][j] + coriolis_velocity[i - 1][j]) * dt;
        }
    }
    return coriolis_drift_distance;
}

// Wrapper for calculate_drift_distances
static PyObject* py_calculate_drift_distances(PyObject* self, PyObject* args) {
    PyObject* py_coriolis_velocity;
    double time;
    int num_steps;

    if (!PyArg_ParseTuple(args, "OdI", &py_coriolis_velocity, &time, &num_steps)) {
        return nullptr;
    }

    // Extract coriolis_velocity from Python list of lists
    std::vector<std::vector<double>> coriolis_velocity;
    if (PyList_Check(py_coriolis_velocity)) {
        for (Py_ssize_t i = 0; i < PyList_Size(py_coriolis_velocity); ++i) {
            PyObject* inner_list = PyList_GetItem(py_coriolis_velocity, i);
            std::vector<double> inner_vector;

            if (PyList_Check(inner_list)) {
                for (Py_ssize_t j = 0; j < PyList_Size(inner_list); ++j) {
                    inner_vector.push_back(PyFloat_AsDouble(PyList_GetItem(inner_list, j)));
                }
            }
            coriolis_velocity.push_back(inner_vector);
        }
    }

    auto drift_distances = calculate_drift_distances(coriolis_velocity, time, num_steps);

    // Build Python list of lists for output
    PyObject* py_result = PyList_New(num_steps);
    for (size_t i = 0; i < drift_distances.size(); ++i) {
        PyObject* py_inner = PyList_New(3);
        for (size_t j = 0; j < 3; ++j) {
            PyList_SetItem(py_inner, j, Py_BuildValue("d", drift_distances[i][j]));
        }
        PyList_SetItem(py_result, i, py_inner);
    }

    return py_result;
}


// Method definitions for the module
static PyMethodDef CoriolisMethods[] = {
    {"radius_at_latitude", py_radius_at_latitude, METH_VARARGS, "Calculate Earth's radius at a given latitude"},
    {"direction_vector", py_direction_vector, METH_VARARGS, "Calculate direction vector between two lat/lon points"},
    {"rotation_matrix", py_rotation_matrix, METH_VARARGS, "Calculate the rotation matrix based on latitude"},
    {"coriolis_acc", py_coriolis_acc, METH_VARARGS, "Calculate Coriolis accelerations"},
    {"calculate_velocities", py_calculate_velocities, METH_VARARGS, "Calculate velocities from Coriolis acceleration"},
    {"calculate_drift_distances", py_calculate_drift_distances, METH_VARARGS, "Calculate drift distances from velocities"},
    {"haversine", py_haversine, METH_VARARGS, "Calculate the great-circle distance using the Haversine formula"},
    {nullptr, nullptr, 0, nullptr} // Sentinel
};

// Module definition
static struct PyModuleDef coriolismodule = {
    PyModuleDef_HEAD_INIT,
    "coriolis_module", // Module name
    nullptr,          // Module documentation
    -1,               // Size of per-interpreter state of the module
    CoriolisMethods   // Methods defined in the module
};

// Module initialization function
PyMODINIT_FUNC PyInit_coriolis_module(void) {
    return PyModule_Create(&coriolismodule);
}
