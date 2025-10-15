# WARP.md

This file provides guidance to WARP (warp.dev) when working with code in this repository.

## Development Commands

### Environment Setup
```bash
# Create and activate virtual environment
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Configure environment
cp .env.example .env
# Edit .env with your Supabase credentials
```

### Running the Application
```bash
# Development server with auto-reload
uvicorn main:app --reload

# Production server
uvicorn main:app --host 0.0.0.0 --port 8000

# Run from main.py directly
python main.py
```

### Testing and Code Quality
```bash
# Run tests
pytest

# Install development tools
pip install black isort flake8

# Format code
black .
isort .

# Lint code
flake8 .
```

## Architecture Overview

### Application Structure
This is a FastAPI-based bioinformatics API with a modular architecture:

- **FastAPI Application**: Main entry point in `main.py`
- **Routers**: API endpoints organized by domain (`app/routers/`)
- **Services**: Business logic layer (`app/services/`)
- **Models**: Pydantic data models for request/response (`app/models/`)
- **Core**: Configuration and database setup (`app/core/`)

### Key Components

#### Database Integration
- Uses **Supabase** as the backend database
- Singleton pattern for database client in `app/core/database.py`
- Configuration managed through environment variables

#### Sequence Analysis Pipeline
The application provides two main analysis domains:

1. **Sequence Operations** (`sequences` router):
   - Basic sequence analysis (length, GC content, composition)
   - Sequence transformations (reverse complement, transcription, translation)
   - Uses Biopython's `Seq` class for molecular biology operations

2. **Annotation Analysis** (`annotations` router):
   - Open Reading Frame (ORF) detection across 6 reading frames
   - Motif pattern matching (TATA box, CAAT box, etc.)
   - GC content analysis with sliding windows

#### Service Layer Architecture
- `SequenceAnalyzer`: Handles sequence analysis using Biopython
- `AnnotationService`: Manages annotation detection and analysis
- Services are stateless and instantiated at router level

### Configuration Management
- Uses Pydantic `BaseSettings` for environment-based configuration
- Required variables: `SUPABASE_URL`, `SUPABASE_KEY`
- Optional: `DEBUG`, `API_TITLE`, `API_VERSION`

### API Design Patterns
- RESTful endpoints with consistent `/api/v1/` prefix
- Pydantic models for request validation and response serialization
- FastAPI automatic documentation at `/docs` and `/redoc`
- Structured error handling with HTTP status codes

### Data Flow
1. Request validation through Pydantic models
2. Business logic processing in service layer
3. Biopython integration for bioinformatics calculations
4. Response serialization and return

## Project-Specific Development Guidelines

### Adding New Endpoints
1. Define request/response models in appropriate `app/models/` file
2. Implement business logic in corresponding service class
3. Create router endpoint with proper error handling
4. Include router in `main.py`

### Working with Sequences
- All sequences are normalized to uppercase for consistency
- Support for DNA, RNA, and protein sequence types
- Use Biopython's `Seq` class for sequence operations

### Error Handling
- Use FastAPI's `HTTPException` for API errors
- Wrap service calls in try-catch blocks at router level
- Return meaningful error messages to clients

### Environment Variables
Never commit actual Supabase credentials. Always use the `.env.example` template and configure your local `.env` file.