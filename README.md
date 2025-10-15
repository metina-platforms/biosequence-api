# BioSequence API

A bioinformatics API for DNA/RNA sequence analysis and annotation identification using FastAPI, Biopython, and Supabase.

## Features

- **Sequence Analysis**: Analyze DNA, RNA, and protein sequences
- **Sequence Manipulation**: Reverse complement, transcription, translation
- **Annotation Detection**: Find Open Reading Frames (ORFs) and motifs
- **GC Content Analysis**: Calculate GC content with sliding window
- **Database Integration**: Store and retrieve data using Supabase
- **RESTful API**: FastAPI-powered endpoints with automatic documentation

## Tech Stack

- **FastAPI**: Modern, fast web framework for building APIs
- **Biopython**: Comprehensive bioinformatics library
- **Supabase**: Backend-as-a-Service for database management
- **Pydantic**: Data validation using Python type annotations
- **Uvicorn**: ASGI web server implementation

## Installation

1. **Clone the repository**
   ```bash
   git clone <repository-url>
   cd biosequence-api
   ```

2. **Create and activate virtual environment**
   ```bash
   python3 -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. **Install dependencies**
   ```bash
   pip install -r requirements.txt
   ```

4. **Environment Configuration**
   ```bash
   cp .env.example .env
   # Edit .env with your Supabase credentials
   ```

## Configuration

Create a `.env` file in the project root:

```env
SUPABASE_URL=your_supabase_url_here
SUPABASE_KEY=your_supabase_anon_key_here
DEBUG=false
API_TITLE=BioSequence API
API_VERSION=1.0.0
```

## Usage

### Starting the Server

```bash
# Development server
uvicorn main:app --reload

# Production server
uvicorn main:app --host 0.0.0.0 --port 8000
```

The API will be available at `http://localhost:8000`

### API Documentation

- **Swagger UI**: `http://localhost:8000/docs`
- **ReDoc**: `http://localhost:8000/redoc`

## API Endpoints

### Sequence Analysis

#### POST `/api/v1/sequences/analyze`
Analyze a biological sequence for basic properties and composition.

**Request Body:**
```json
{
  "sequence": "ATGCGTACGTAGCTA",
  "sequence_type": "DNA"
}
```

**Response:**
```json
{
  "sequence": "ATGCGTACGTAGCTA",
  "sequence_type": "DNA",
  "length": 15,
  "gc_content": 53.33,
  "molecular_weight": 4563.98,
  "composition": {
    "A": 4,
    "T": 4,
    "G": 4,
    "C": 3
  }
}
```

#### POST `/api/v1/sequences/reverse-complement`
Get the reverse complement of a DNA sequence.

#### POST `/api/v1/sequences/transcribe`
Transcribe DNA sequence to RNA.

#### POST `/api/v1/sequences/translate`
Translate DNA/RNA sequence to protein.

### Annotation Analysis

#### POST `/api/v1/annotations/find-orfs`
Find Open Reading Frames in a DNA sequence.

**Request Body:**
```json
{
  "sequence": "ATGAAATTTAAATAG",
  "min_length": 12
}
```

#### POST `/api/v1/annotations/find-motifs`
Find common sequence motifs (TATA box, CAAT box, etc.).

#### POST `/api/v1/annotations/gc-content`
Calculate GC content using sliding window analysis.

## Project Structure

```
biosequence-api/
├── app/
│   ├── core/
│   │   ├── __init__.py
│   │   ├── config.py          # Application configuration
│   │   └── database.py        # Supabase client setup
│   ├── models/
│   │   ├── __init__.py
│   │   ├── sequence.py        # Sequence data models
│   │   └── annotation.py      # Annotation data models
│   ├── routers/
│   │   ├── __init__.py
│   │   ├── sequences.py       # Sequence analysis endpoints
│   │   └── annotations.py     # Annotation analysis endpoints
│   ├── services/
│   │   ├── __init__.py
│   │   ├── sequence_service.py    # Sequence analysis logic
│   │   └── annotation_service.py  # Annotation analysis logic
│   └── __init__.py
├── tests/
├── venv/
├── .env.example
├── .gitignore
├── main.py                    # Application entry point
├── requirements.txt
└── README.md
```

## Development

### Running Tests

```bash
pytest
```

### Code Formatting

```bash
# Install development dependencies
pip install black isort flake8

# Format code
black .
isort .

# Check linting
flake8 .
```

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add some amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- [Biopython](https://biopython.org/) for bioinformatics tools
- [FastAPI](https://fastapi.tiangolo.com/) for the web framework
- [Supabase](https://supabase.com/) for the backend infrastructure