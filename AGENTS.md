# AGENTS.md

This file provides essential instructions, build commands, and coding standards for agentic coding tools operating in this repository.

## 🛠 Build & Test Commands

### Building
The project uses **CMake**. All builds should be performed in a dedicated `build` directory.

- **General Build (All Platforms):**
  ```bash
  mkdir build
  cd build
  cmake ..
  cmake --build . --config Release
  ```

- **Visual Studio (Windows):**
  Run `create_sln.bat` in the root, then open `build/box2d.sln`.

- **Linux:**
  Run `build.sh` from the root.

- **macOS (Xcode):**
  ```bash
  mkdir build
  cd build
  cmake -G Xcode ..
  ```

### Testing
Tests are enabled by default (`BOX2D_UNIT_TESTS=ON`).

- **Run All Tests:**
  ```bash
  cd build
  cmake --build . --target test
  ```

### Linting
The project can treat warnings as errors. When contributing, ensure your code passes with this enabled:
```bash
cmake .. -DBOX2D_COMPILE_WARNING_AS_ERROR=ON
```

## 📜 Code Style & Guidelines

### Language Standards
- **Core Library:** C17.
- **Samples & Tests:** C++20.

### Naming Conventions
- **Files:** `snake_case.c` and `snake_case.h`.
- **Functions/Variables:** standard C naming (lowercase with underscores).

### Coding Patterns
- **Data-Oriented Design:** Focus on memory layout and cache efficiency.
- **SIMD:** The library utilizes SSE2 and NEON. Do not break SIMD-friendly structures unless explicitly intended.
- **Memory Management:** Uses arena allocators (`arena_allocator.h`) for efficiency. Use them where appropriate to avoid fragmentation.
- **Error Handling:** Heavy use of validation logic. Ensure `BOX2D_VALIDATE` is respected during development.

### Architecture
- **Threading:** Uses `enkiTS` for task-based parallelism in samples and tests.
- **Optimization:** Code is highly optimized for performance. Avoid introducing overhead that disrupts throughput or latency.
