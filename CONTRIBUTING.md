# Contributing to SpatialRegret

Thank you for your interest in contributing to the SpatialRegret project! This document provides guidelines for contributing to the codebase.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [How Can I Contribute?](#how-can-i-contribute)
- [Development Setup](#development-setup)
- [Coding Standards](#coding-standards)
- [Submitting Changes](#submitting-changes)
- [Reporting Bugs](#reporting-bugs)
- [Suggesting Enhancements](#suggesting-enhancements)

## Code of Conduct

This project adheres to a code of conduct that we expect all contributors to follow:

- Be respectful and inclusive
- Welcome newcomers and help them get started
- Focus on constructive feedback
- Accept responsibility and apologize for mistakes

## How Can I Contribute?

### Types of Contributions

We welcome several types of contributions:

1. **Bug Reports**: Found a bug? Let us know!
2. **Feature Requests**: Have an idea for improvement? Share it!
3. **Code Contributions**: Fix bugs, add features, or improve performance
4. **Documentation**: Improve README, examples, or function documentation
5. **Testing**: Add test cases or validate existing functionality

### Good First Issues

If you're new to the project, look for issues labeled:
- `good first issue`
- `documentation`
- `help wanted`

## Development Setup

1. **Fork the repository** on GitHub

2. **Clone your fork**:
   ```bash
   git clone https://github.com/YOUR-USERNAME/SpatialRegret.git
   cd SpatialRegret
   ```

3. **Add upstream remote**:
   ```bash
   git remote add upstream https://github.com/DecodEPFL/SpatialRegret.git
   ```

4. **Install dependencies** (see [REQUIREMENTS.md](REQUIREMENTS.md))

5. **Create a branch** for your changes:
   ```bash
   git checkout -b feature/your-feature-name
   ```

## Coding Standards

### MATLAB Code Style

Please follow these conventions:

#### Naming Conventions
- **Functions**: `lowercase_with_underscores.m`
- **Variables**: `camelCase` or `lowercase_with_underscores`
- **Constants**: `UPPERCASE_WITH_UNDERSCORES`
- **Classes**: `PascalCase`

#### Code Structure
```matlab
function [outpuP11, outpuP12] = function_name(inpuP11, inpuP12, options)
% FUNCTION_NAME - Brief description (one line)
%
% Detailed description of what the function does, explaining the 
% algorithm or methodology if relevant.
%
% Syntax:
%   [outpuP11, outpuP12] = function_name(inpuP11, inpuP12)
%   [outpuP11, outpuP12] = function_name(inpuP11, inpuP12, options)
%
% Inputs:
%   inpuP11  - Description of inpuP11 (type, dimensions, units)
%   inpuP12  - Description of inpuP12
%   options - (optional) Structure with optional parameters
%             .field1 - Description
%             .field2 - Description
%
% Outputs:
%   outpuP11 - Description of outpuP11
%   outpuP12 - Description of outpuP12
%
% Example:
%   delays = eye(5);
%   [K, Q, cost] = function_name(sys, delays);
%
% See also: RELATED_FUNCTION1, RELATED_FUNCTION2

% Input validation
if nargin < 2
    error('At least 2 inputs required');
end

if nargin < 3
    options = struct();  % Default options
end

% Main code here
...

end
```

#### Comments
- Add comments for complex logic
- Use `%` for single-line comments
- Use `%%` for section breaks
- Document all function inputs/outputs

#### Code Formatting
- **Indentation**: 4 spaces (no tabs)
- **Line length**: Try to keep lines under 80 characters
- **Whitespace**: Use blank lines to separate logical sections

### Performance Considerations

- Preallocate arrays when possible
- Avoid nested loops for large data
- Use vectorized operations when appropriate
- Clear large variables when no longer needed

### Testing

Before submitting:

1. **Run your code** on a small example
2. **Verify results** are correct
3. **Check for warnings** 
4. **Test edge cases** (empty inputs, large systems, etc.)

Example test pattern:
```matlab
% Test with small system
n_agents = 3;
... % setup
% Verify output dimensions
assert(size(K, 1) == expected_rows, 'Incorrect output size');
```

## Submitting Changes

### Pull Request Process

1. **Update your branch** with latest upstream:
   ```bash
   git fetch upstream
   git rebase upstream/main
   ```

2. **Test your changes** thoroughly

3. **Commit with clear messages**:
   ```bash
   git commit -m "Add feature: brief description
   
   Detailed explanation of what changed and why.
   Fixes #issue-number"
   ```

4. **Push to your fork**:
   ```bash
   git push origin feature/your-feature-name
   ```

5. **Create a Pull Request** on GitHub

### Pull Request Guidelines

Your PR should:

- Have a **clear title** describing the change
- Include a **description** of what changed and why
- Reference any related **issues** (e.g., "Fixes #123")
- Include **examples** or **tests** if applicable
- Update **documentation** if needed

### PR Checklist

- [ ] Code follows the style guidelines
- [ ] Comments and documentation are updated
- [ ] Changes have been tested
- [ ] No new warnings introduced
- [ ] Examples still run correctly
- [ ] README updated if needed

## Reporting Bugs

### Before Submitting a Bug Report

- Check existing issues to avoid duplicates
- Try to reproduce with the latest version
- Collect relevant information (MATLAB version, OS, error messages)

### Bug Report Template

```markdown
**Description**
A clear description of the bug.

**To Reproduce**
Steps to reproduce the behavior:
1. Run script '...'
2. With parameters '...'
3. See error

**Expected Behavior**
What you expected to happen.

**Actual Behavior**
What actually happened, including error messages.

**Environment**
- MATLAB Version: [e.g., R2021a]
- OS: [e.g., macOS 12.0]
- Toolboxes: [e.g., Control System Toolbox v10.9]
- Solvers: [e.g., Gurobi 9.5]

**Additional Context**
Any other relevant information.
```

## Suggesting Enhancements

We welcome feature suggestions! Please:

1. **Check existing issues** for similar suggestions
2. **Describe the use case** and problem it solves
3. **Explain your proposed solution**
4. **Consider alternatives** you've thought about

### Enhancement Template

```markdown
**Feature Description**
Clear description of the proposed feature.

**Motivation**
Why is this feature needed? What problem does it solve?

**Proposed Solution**
How would you implement this?

**Alternatives Considered**
Other approaches you've thought about.

**Additional Context**
Any other relevant information.
```

## Questions?

If you have questions about contributing:

- Open an issue with the `question` label
- Contact the maintainers via the GitHub repository

## Attribution

By contributing, you agree that your contributions will be licensed under the same license as the project (CC by 4.0).

## Thank You!

Your contributions make this project better for everyone. We appreciate your time and effort! 🎉
