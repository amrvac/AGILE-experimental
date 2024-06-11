#!/usr/bin/env python3
from kernel_tuner import tune_kernel, run_kernel
from kernel_tuner.utils.directives import (
    Code,
    OpenACC,
    Fortran,
    extract_directive_signature,
    extract_directive_code,
    extract_initialization_code,
    extract_preprocessor,
    generate_directive_function,
    extract_directive_data,
    allocate_signature_memory,
)


if __name__ == '__main__':
    with open('laplace_kernel_tuner.f90') as f:
        code = f.read()

    name = "update_solution"
    sizes = {'n_vars': 2,
             'grid_size': 1024,
             'block_size': 16,
             'n_blocks': 64}

    # Extract tunable directive
    app = Code(OpenACC(), Fortran())
    preprocessor = extract_preprocessor(code)

    # add user dims to preprocessor
    for k, v in sizes.items():
        preprocessor.append(f"#define {k} {v}")

    signature = extract_directive_signature(code, app)
    body = extract_directive_code(code, app)
    # Allocate memory on the host
    data = extract_directive_data(code, app)
    init = extract_initialization_code(code, app)
    args = allocate_signature_memory(data[name], preprocessor) #, user_dimensions=sizes)
    # Generate kernel string
    kernel_string = generate_directive_function(
        preprocessor, signature[name], body[name], app, data=data[name], initialization=init
    )
    # print("preprocessor:")
    # print(preprocessor)
    # print()
    # print("signature:")
    # print(signature[name])
    # print()
    # print("body:")
    # print(body[name])
    # print()
    print("kernel string:")
    print(kernel_string)
    print()

    tune_params = {'collapse_factor': [1, 2, 3],
                   'vlength': [2**i for i in range(0, 11, 2)]}

    defines = {k: k for k in tune_params}

    compiler_options = ["-acc", "-fast", "-Mpreprocess"]
    tune_kernel(name, kernel_string, 0, args, tune_params, compiler_options=compiler_options, compiler="nvfortran", defines=defines)
