from typing import Optional
from supabase import create_client, Client
from app.core.config import settings


class SupabaseClient:
    _client: Optional[Client] = None
    
    @classmethod
    def get_client(cls) -> Client:
        if cls._client is None:
            if not settings.supabase_url or not settings.supabase_key:
                raise ValueError("Supabase URL and key must be configured")
            cls._client = create_client(settings.supabase_url, settings.supabase_key)
        return cls._client
    
    @classmethod
    def close(cls):
        if cls._client:
            cls._client = None


# Convenience function to get the client
def get_supabase_client() -> Client:
    return SupabaseClient.get_client()