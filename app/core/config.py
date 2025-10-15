import os
from typing import Optional
from pydantic_settings import BaseSettings


class Settings(BaseSettings):
    # Supabase configuration
    supabase_url: Optional[str] = None
    supabase_key: Optional[str] = None
    
    # Database configuration
    database_url: Optional[str] = None
    
    # API configuration
    api_title: str = "BioSequence API"
    api_version: str = "1.0.0"
    debug: bool = False
    
    class Config:
        env_file = ".env"


settings = Settings()