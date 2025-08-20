#!/usr/bin/env python3
"""
Repository Policy Self-Checker

Validates repository structure against .github/REPOSITORY_POLICY.md
ensuring compliance with established standards.
"""

import os
import sys
from pathlib import Path
from typing import List, Dict, Tuple


class PolicyChecker:
    """Checks repository compliance with REPOSITORY_POLICY.md standards."""
    
    def __init__(self, repo_root: str = "."):
        self.repo_root = Path(repo_root)
        self.issues: List[str] = []
        self.warnings: List[str] = []
        
    def check_required_directories(self) -> bool:
        """Check that required directories exist."""
        required_dirs = [
            ".github",
            ".grok", 
            "applications",
            "data",
            "docs",
            "notebooks", 
            "proof_pack",
            "static",
            "templates",
            "tests",
            "wave-ribonn"
        ]
        
        success = True
        for dirname in required_dirs:
            dir_path = self.repo_root / dirname
            if not dir_path.exists():
                self.issues.append(f"Missing required directory: {dirname}")
                success = False
            elif not dir_path.is_dir():
                self.issues.append(f"Required path is not a directory: {dirname}")
                success = False
                
        return success
    
    def check_required_files(self) -> bool:
        """Check that policy-required files exist."""
        required_files = [
            ".github/REPOSITORY_POLICY.md",
            ".github/copilot-instructions.yml",
            "requirements.txt",
            "README.md",
            "run_tests.py"
        ]
        
        success = True
        for filepath in required_files:
            file_path = self.repo_root / filepath
            if not file_path.exists():
                self.issues.append(f"Missing required file: {filepath}")
                success = False
            elif not file_path.is_file():
                self.issues.append(f"Required path is not a file: {filepath}")
                success = False
                
        return success
    
    def check_naming_conventions(self) -> bool:
        """Check file naming conventions."""
        success = True
        
        # Check Python files are snake_case
        for py_file in self.repo_root.glob("*.py"):
            if not self._is_snake_case(py_file.stem):
                self.warnings.append(f"Python file not snake_case: {py_file.name}")
        
        # Check major docs are UPPERCASE
        major_docs = ["README.md", "LICENSE"]
        for doc in self.repo_root.glob("*.md"):
            if doc.name.startswith("README") or doc.name.startswith("LICENSE"):
                continue
            if "_" in doc.stem and not doc.stem.isupper():
                # Allow mixed case for issue resolution docs
                if "ISSUE_" not in doc.name and "SUMMARY" not in doc.name:
                    self.warnings.append(f"Major doc not UPPERCASE: {doc.name}")
                    
        return success
    
    def check_core_modules(self) -> bool:
        """Check that core scientific modules exist at root level."""
        expected_modules = [
            "z_framework.py",
            "topological_analysis.py", 
            "invariant_features.py"
        ]
        
        success = True
        for module in expected_modules:
            module_path = self.repo_root / module
            if not module_path.exists():
                self.issues.append(f"Missing core module: {module}")
                success = False
                
        return success
    
    def check_test_structure(self) -> bool:
        """Check test organization (dual location allowance)."""
        success = True
        
        # Check for root-level test files
        root_tests = list(self.repo_root.glob("test_*.py"))
        tests_dir_tests = list((self.repo_root / "tests").glob("test_*.py"))
        
        if not root_tests and not tests_dir_tests:
            self.issues.append("No test files found in root or tests/ directory")
            success = False
        
        # Check run_tests.py orchestrator exists
        if not (self.repo_root / "run_tests.py").exists():
            self.issues.append("Missing test orchestrator: run_tests.py")
            success = False
            
        return success
    
    def check_dependencies(self) -> bool:
        """Check requirements.txt has exact version pinning."""
        success = True
        req_file = self.repo_root / "requirements.txt"
        
        if not req_file.exists():
            return False
            
        try:
            with open(req_file, 'r') as f:
                lines = f.readlines()
                
            for line in lines:
                line = line.strip()
                if line and not line.startswith('#'):
                    if '==' not in line:
                        self.warnings.append(f"Dependency not pinned: {line}")
                        
        except Exception as e:
            self.issues.append(f"Could not read requirements.txt: {e}")
            success = False
            
        return success
    
    def _is_snake_case(self, name: str) -> bool:
        """Check if string is snake_case."""
        if not name:
            return False
        return name.lower() == name and '_' in name and ' ' not in name
    
    def run_all_checks(self) -> bool:
        """Run all policy checks."""
        print("🔍 Running repository policy checks...")
        print("=" * 50)
        
        checks = [
            ("Required directories", self.check_required_directories),
            ("Required files", self.check_required_files), 
            ("Naming conventions", self.check_naming_conventions),
            ("Core modules", self.check_core_modules),
            ("Test structure", self.check_test_structure),
            ("Dependencies", self.check_dependencies)
        ]
        
        all_passed = True
        for check_name, check_func in checks:
            try:
                result = check_func()
                status = "✓ PASS" if result else "❌ FAIL"
                print(f"{check_name}: {status}")
                if not result:
                    all_passed = False
            except Exception as e:
                print(f"{check_name}: ❌ ERROR - {e}")
                all_passed = False
        
        # Report issues and warnings
        if self.issues:
            print(f"\n❌ ISSUES FOUND ({len(self.issues)}):")
            for issue in self.issues:
                print(f"  • {issue}")
                
        if self.warnings:
            print(f"\n⚠️  WARNINGS ({len(self.warnings)}):")
            for warning in self.warnings:
                print(f"  • {warning}")
        
        print("=" * 50)
        if all_passed and not self.issues:
            print("🎉 All policy checks PASSED!")
            return True
        else:
            print("💥 Some policy checks FAILED!")
            return False


def main():
    """Main entry point."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Check repository policy compliance")
    parser.add_argument("--repo-root", default=".", 
                      help="Repository root directory (default: current directory)")
    
    args = parser.parse_args()
    
    checker = PolicyChecker(args.repo_root)
    success = checker.run_all_checks()
    
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()